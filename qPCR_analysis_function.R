library(ggplot2)
library(viridis)

#### FUNCTIONS ####

qpcr_analysis <- function( # MAIN FUNCTION
  ct_data, # Data frame, expected columns: names, cq, efficiency, primers, included
  type = "rtqpcr", # Analysis type, rtqpcr or qchip
  design, # Design matrix, at least one expected column "sample_name" and additional columns for each exp condition, for chip at least one has to be named IP
  calibsample, # String, name of the reference sample for dct calculation, will be 1 for all genes and the rest of the samples will be relative to it
  hkg, # Character vector with the "housekeeping" genes for normalization
  exp_name, # String, experiment name, will be used for naming saved files
  fix_names = FALSE, # Logical, fix names from chainy output if TRUE
  exclude = FALSE, # Logical, exclude samples with included = FALSE
  save_csv = TRUE # Logical save a csv file with the normalized data
) {
  
  # Fix names if ncessary (chainy output), otherwise keep as is
  
  if (fix_names) {
    ct_data$sample <- sapply(as.character(ct_data$names), fix_names, USE.NAMES = F)
  } else {
    ct_data$sample <- ct_data$names
  }
  
  # Exclude samples if necessary
  
  if (exclude) {
    ct_data <- ct_data[ct_data$included == TRUE,]
  }
  
  # Average ct of technical replicates
  
  ct_avg <- ct_avg_fun(ct_data)
  
  # Calculate median efficieny and average ct per primer pair
  
  eff <- aggregate(ct_avg$effaverage, by = ct_avg["primers"], FUN = median, na.rm = T)
  names(eff)[2] <- "efficiency"
  print(eff)
  ct <- aggregate(ct_data$cq, by = ct_data["primers"], FUN = mean, na.rm = T)
  names(ct)[2] <- "cq"
  save(eff, file =  paste0(exp_name, "_eff.RData"))
  write.csv(eff, paste0(exp_name, "_eff.csv"))
  save(ct, file = paste0(exp_name, "_ct.RData"))
  write.csv(ct, paste0(exp_name, "_ct.csv"))
  
  
  if (type == "rtqpcr") {
    
    # Calculate dct for the data
    
    dct_data <- dct_apply(ct_avg, calibsample, eff)
    
    # Normalization factor calculation
    
    norm_factor <- norm_factor_fun(dct_data, hkg, exp_name)
    
    # Normalize dct values
    
    norm_dct <- apply(dct_data, 1, function(x) as.numeric(x["dct"]) / norm_factor[as.character(x["sample"])])
    norm_data <- cbind(dct_data, "norm_dct" = norm_dct)
    write.csv(norm_data, file = paste0(exp_name, "_norm_data.csv"))
    
    # Create QC plots
    ct_sd_plot(norm_data, exp_name = exp_name)
    ct_avg_plot(norm_data, exp_name = exp_name)
    eff_avg_plot(norm_data, exp_name = exp_name)
    
    return(norm_data)
  } 
}

fix_names <- function(x) { # fix names from chainy output (remove plate well coordinates)
  name <- as.character(x)
  index <- unlist(gregexpr("..", x, fixed = T)) + 2
  indexstar <- unlist(gregexpr("**", x, fixed = T))
  if (indexstar == -1) {
    fixed_name <- substr(name, index, nchar(x))
  } else {
    fixed_name <- substr(name, index, indexstar[2] - 1)
  }
  fixed_name
}

ct_avg_fun <- function(ct_data) {
  ct_avg <- data.frame()
  for (sample in unique(ct_data$sample)) {
    sample_data <- ct_data[ct_data$sample == sample,]
    sample_name <- sample_data[,"sample"]
    avg_ct <- aggregate(sample_data$cq, by = sample_data["primers"], FUN = mean, na.rm = F)
    names(avg_ct) <- c("primers", "ctaverage")
    avg_eff <- aggregate(sample_data$efficiency, by = sample_data["primers"], FUN = mean, na.rm = F)
    names(avg_eff) <- c("primers", "effaverage")
    ct_sd <- aggregate(sample_data$cq, by = sample_data["primers"], FUN = sd, na.rm = F)
    names(ct_sd) <- c("primers", "ct_sd")
    avg_sample_data <- data.frame("sample" = rep(sample_data[1, "sample"], nrow(avg_ct)),
                                  "primers" = avg_ct$primers, "effaverage" = avg_eff$effaverage, "ctaverage" = avg_ct$ctaverage, "ct_sd" = ct_sd$ct_sd)
    ct_avg <- rbind(ct_avg, avg_sample_data)
  }
  ct_avg
}

norm_factor_fun <- function(dct_data, hkg, exp_name) { # calculate normalization factor for each sample and generate a scatterplot of housekeeping gene expression
  hkg_data <- dct_data[grep(paste(hkg, collapse = "|"), dct_data$primers),]
  hkg_centered <- aggregate(hkg_data$dct, by = hkg_data["primers"], FUN = function(x) x)
  hkg_centered <- t(hkg_centered[,2:ncol(hkg_centered)])
  colnames(hkg_centered) <- hkg
  rownames(hkg_centered) <- unique(hkg_data$sample)
  norm_factor <- apply(hkg_centered, 1, mean, na.rm = T)
  
  # Scatter plot of HK gene expression values
  hkg_centered <- as.data.frame(hkg_centered)
  hkg_plot <- ggplot(hkg_centered, aes(x = get(hkg[1]), y = get(hkg[2]), label = rownames(hkg_centered))) +
    coord_cartesian(xlim = c(min(hkg_centered), max(hkg_centered)),
                    ylim = c(min(hkg_centered), max(hkg_centered))) +
    geom_abline(slope = 1, intercept = 0, color = "darkblue") +
    geom_point() +
    geom_text(hjust = 0, nudge_x = 0.05) +
    labs(x = names(hkg_centered)[1], 
         y = names(hkg_centered)[2])
  
  ggsave(filename = paste0(exp_name, "_hkg_scatter.pdf"), 
         plot = hkg_plot, 
         device = cairo_pdf) 
  
  norm_factor
}

dct <- function(x, calib, eff.matrix) { # takes a row of data for a single sample, the ct values
  # of the calibrator sample and the efficiency of each primer pair, computes the dct value
  ct <- as.numeric(x["ctaverage"])
  sample <- x["sample"]
  primers <- x["primers"]
  eff.primers <- eff.matrix[eff.matrix$primers == primers, "efficiency"]
  dct_value <- eff.primers^(calib[primers] - ct)
  dct_value
}

dct_apply <- function (data, calibsample, eff) { # helper function to apply the dct function to each row
  ctcalib <- data[data$sample == calibsample, "ctaverage"]
  names(ctcalib) <- data[data$sample == calibsample, "primers"]
  dct_values <- apply(data, 1, dct, calib = ctcalib, eff.matrix = eff)
  new_data <- cbind(data, "dct" = dct_values)
  new_data
}

dct_apply_chip <- function (data, ...) { # helper function to apply the dct function to each row
  ctcalib <- data[data$input, "ctaverage"]
  names(ctcalib) <- unique(data$primers)
  dct_values <- apply(data, 1, dct, calib = ctcalib, eff.matrix = eff)
  data$dct <- dct_values
  data
}

dct_group_chip <- function(dct_data, ...) {
  dct_data <- NULL
  for (i in 1:length(unique(ct_avg$group))) {
    group <- unique(ct_avg$group)[i]
    data <- ct_avg[ct_avg$group == group,]
    new_data <- dct_apply_chip(data, eff)
    if (is.null(dct_data)) {
      dct_data <- new_data
    } else {
      dct_data <- rbind(dct_data, new_data)
    }
  }
  dct_data
}

read_chainy_dir <- function(directory = "./chainy_output/") {
  file_list <- list.files(directory)
  print(paste0("Files found: ", file_list))
  chainy_output <- data.frame()
  for (file in file_list) {
    temp_data <- read.delim(paste(directory, file, sep = "/"))
    index <- unlist(gregexpr("_cq_eff", file, fixed = T)) - 1
    target <- rep(substr(file, 0, index), nrow(temp_data))
    target <- toupper(target) # comment if you don't want to get all primer names in uppercase
    temp_data <- cbind(temp_data, "primers" = target)
    chainy_output <- rbind(chainy_output, temp_data)
  }
  names(chainy_output) <- c("names", "cq", "efficiency", "primers")
  chainy_output
}

batch_correction_chip <- function(data, reference, value) {
    # Expects the following columns in the data:
    # "group" with your unique experimental groups
    # "biorep" with the unique experimental replicates (it is not necessary numeric)
  if (!"batch_cor" %in% names(data)) {
    data <- cbind(data, "batch_cor" = rep(0, length.out = nrow(data)))
  }
  for (i in unique(data$primers)) {
    for (j in unique(data$biorep)) {
      target_data <- data[data$primers == i & data$biorep == j, value]
      data[data$primers == i & data$biorep == j, "batch_cor"] <- target_data / mean(target_data, na.rm = T)
    }
    reference_avg <- mean(data[data$primers == i & data$group == reference, "batch_cor"], na.rm = T)
    data[data$primers == i, "batch_cor"] <- data[data$primers == i, "batch_cor"] / reference_avg
  }
  data
}

ct_sd_plot <- function(norm_data, max_sd = sd(c(0, 0.5)), max_ct = 35, exp_name = "exp_name") {
  ctsd <- ggplot(norm_data, aes(x = ctaverage, y = ct_sd, color = primers, fill = primers)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = max_sd, color = "red") +
    geom_vline(xintercept = max_ct, color = "red") +
    scale_color_viridis(option = "turbo", discrete = TRUE)
  ggsave(filename = paste0(exp_name, "_ct_sd_plot.pdf"), 
         plot = ctsd, 
         device = cairo_pdf) 
}

ct_avg_plot <- function(norm_data, max_ct = 35, exp_name = "exp_name") {
  ctavg <- ggplot(norm_data, aes(x = primers, y = ctaverage, color = primers, group = primers, size = ct_sd)) +
    geom_jitter(width = 0.2, alpha = 0.75) +
    geom_hline(yintercept = max_ct, color = "red") +
    scale_color_viridis(option = "turbo", discrete = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0(exp_name, "_ct_avg_plot.pdf"), 
         plot = ctavg, 
         device = cairo_pdf) 
}

eff_avg_plot <- function(norm_data, effrange = c(1.8, 2), exp_name = "exp_name") {
  effavg <- ggplot(norm_data, aes(x = primers, y = effaverage, color = ctaverage, group = primers)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = effrange, color = "red") +
    scale_color_viridis(option = "turbo")
  ggsave(filename = paste0(exp_name, "_eff_avg_plot.pdf"), 
         plot = effavg, 
         device = cairo_pdf) 
}
