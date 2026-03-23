####################################################################################################
#                 DRY_NAM: Extraction of drought facets in a bimodal context
#
#                 more information on: !!!
####################################################################################################

# Loading necessary packages
library(parallel)
library(dplyr)
library('phenofit')
library(tidyr)
library(lubridate)
#library(RamanMP)
library(terra)

#setwd("C:/Users/UM/DRY_NAM/SWC")
setwd("C:/Users/bouss/Desktop/PhD/Data_Monsoon")

# Define the output directory
output_dir <- "C:/Users/bouss/Desktop/PhD/Data_Monsoon/"
  
# Load the files to process
#files_to_process <-  c("")

#### 1. Initializing the data frame ####

date <- seq(as.Date("1960-01-01",format="%Y-%m-%d"),as.Date("2024-12-31",format="%Y-%m-%d"),by="day")

years <- format(date, "%Y")
number_of_years <- length(unique(years))
seq_years <- seq(min(as.numeric(years)), max(as.numeric(years)))

DF <- as.data.frame(matrix(NA, nrow = number_of_years,ncol = 45))
colnames(DF) <- c("Year", "First.Drying.Rate", "First.Wetting.Rate",
                  "First.Intercept.drying", "First.Intercept.wetting", "First.Peak.Value",
                  "First.Peak.DOY", "First.Low.D.Onset",
                  "First.Mod.D.Onset", "First.Extreme.D.Onset", "First.Low.D.Offset", "First.Mod.D.Offset",
                  "First.Extreme.D.Offset", "First.Low.D.Duration", "First.Mod.D.Duration",
                  "First.Extreme.D.Duration", "First.Low.D.S", "First.Mod.D.S", "First.Extreme.D.S","Second.Drying.Rate", "Second.Wetting.Rate",
                  "Second.Intercept.drying", "Second.Intercept.wetting", "Second.Peak.Value",
                  "Second.Peak.DOY","Second.Low.D.Onset",
                  "Second.Mod.D.Onset", "Second.Extreme.D.Onset", "Second.Low.D.Offset", "Second.Mod.D.Offset",
                  "Second.Extreme.D.Offset", "Second.Low.D.Duration", "Second.Mod.D.Duration",
                  "Second.Extreme.D.Duration", "Second.Low.D.S", "Second.Mod.D.S", "Second.Extreme.D.S", 
                  "Wetting.Phase.Duration", "Dry.Monsoon.Duration", "Mod.Monsoon.Duration","Wet.Monsoon.Duration", "Peak.Monsoon", "Peak.Monsoon.DOY", 
                  "Cumulative.Year", "Unimodal")


DF$Year <- seq_years

# Drought index parameters
max_DI <- 100 # Drought index limit value
low_onset_threshold <- low_offset_threshold <- max_DI * 0.25
mod_onset_threshold <- mod_offset_threshold <- max_DI * 0.50
high_onset_threshold <- high_offset_threshold <- max_DI * 0.75

#### 2. Load the functions ####

# Peak detection function derived from the "RamanMP" package (Nava V., Frezzotti M. L., Leoni B. (2021), Raman spectroscopy for the analysis of microplastics in aquatic systems. Applied Spectroscopy, 75(11), 1341-1357.). 
# The plots are removed to accelerate the data processing.
peakfunction <- function (spectrum, threshold = 8, m = 20, max.peak = 3) {
  spectrum <- na.omit(spectrum)
  shape <- diff(sign(diff(spectrum[, 2], na.pad = FALSE)))
  colnames(spectrum) <- c("Var1", "Var2")
  spectrum_sin <- spectrum[, 2]
  pks <- sapply(which(shape < 0), FUN = function(i) {
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(spectrum_sin), w, length(spectrum_sin))
    if (all(spectrum_sin[c(z:i, (i + 2):w)] <= spectrum_sin[i + 
                                                            1])) 
      return(i + 1)
    else return(numeric(0))
  })
  peak_index <- unlist(pks)
  peak_inten <- spectrum_sin[peak_index]
  spectra_new <- subset(spectrum, spectrum[, 2] > threshold)
  value_x <- spectra_new[spectra_new[, 2] %in% peak_inten, 
  ]
  min_yaxis <- min(spectrum[, 2])
  max_yaxis <- max(spectrum[, 2])
  diff_yaxis <- max_yaxis - min_yaxis
  gap <- diff_yaxis/2
  nudge_value <- diff_yaxis/2
  Var1 <- Var2 <- NULL
  if (max.peak == 0) {
    value <- value_x[, 1]
    # plot_peaks <- ggplot2::ggplot(spectrum, ggplot2::aes(x = Var1, 
    #                                                      y = Var2, label = round(Var1, digits = 0))) + ggplot2::geom_path(col = "gray5") + 
    #   ggrepel::geom_text_repel(data = spectrum %>% dplyr::filter(spectrum[, 
    #                                                                       1] %in% value), nudge_y = nudge_value, segment.size = 0.3, 
    #                            segment.color = "grey50", direction = "both", 
    #                            angle = 90) + ggplot2::ylim(min_yaxis, (max_yaxis + 
    #                                                                      gap)) + ggplot2::labs(x = expression(Raman ~ shift ~ 
    #                                                                                                             (cm^-1)), y = "Intensity (a.u.)") + ggplot2::theme_minimal() + 
    #   ggplot2::theme(axis.title.x = ggplot2::element_text(size = "16"), 
    #                  axis.title.y = ggplot2::element_text(size = "16"), 
    #                  axis.text.x = ggplot2::element_text(size = "14"), 
    #                  axis.text.y = ggplot2::element_blank())
    # print(plot_peaks)
    colnames(value_x) <- c("Wavelength", "Intensity")
    return(value_x)
  }
  else {
    new_order <- value_x[order(value_x[, 2], decreasing = TRUE), 
    ]
    new_selection <- new_order[1:max.peak, ]
    value2 <- new_selection[, 1]
    # plot_peaks2 <- ggplot2::ggplot(spectrum, ggplot2::aes(x = Var1, 
    #                                                       y = Var2, label = round(Var1, digits = 0))) + ggplot2::geom_line(col = "gray5") + 
    #   ggrepel::geom_text_repel(data = spectrum %>% dplyr::filter(spectrum[, 
    #                                                                       1] %in% value2), nudge_y = nudge_value, segment.size = 0.3, 
    #                            segment.color = "grey50", direction = "both", 
    #                            angle = 90) + ggplot2::ylim(min_yaxis, (max_yaxis + 
    #                                                                      gap)) + ggplot2::labs(x = expression(Raman ~ shift ~ 
    #                                                                                                             (cm^-1)), y = "Intensity (a.u.)") + ggplot2::theme_minimal() + 
    #   ggplot2::theme(axis.title.x = ggplot2::element_text(size = "16"), 
    #                  axis.title.y = ggplot2::element_text(size = "16"), 
    #                  axis.text.x = ggplot2::element_text(size = "14"), 
    #                  axis.text.y = ggplot2::element_blank())
    # print(plot_peaks2)
    new_final <- new_selection[order(new_selection[, 1], 
                                     decreasing = FALSE), ]
    colnames(new_final) <- c("Wavelength", "Intensity")
    return(new_final)
  }
}

# Main DRY_NAM function for drought facets extraction
process_file <- function(file, date, seq_years, DF, low_onset_threshold, low_offset_threshold,
                         mod_onset_threshold,mod_offset_threshold,
                         high_onset_threshold, high_offset_threshold, output_dir) {
  # Load necessary packages
  library(terra)
  library(lubridate)
  library(dplyr)
  library('phenofit')
  library(tidyr)
  
  # Read the csv file
  res.model <- read.csv(file, sep = ";", header = TRUE)
  res.model$date <- as.Date(date)
  res.model$DOY <- yday(res.model$date)
  res.model$year <- year(res.model$date)
  
  # Field capacity value necessary to convert SWC model-based values into the SOLDI index
  FC_tot <-  145.9156
  
  # Converting soil water content into the SOLDI index
  colnames(res.model)[colnames(res.model) == "valeurs"] <- "SOLDI"
  res.model$SOLDI <- - res.model$SOLDI + FC_tot
  res.model$SOLDI[res.model$SOLDI < 0] <- 0
  
  # Retrieve the pixel coordinates
  coord_x_y <- gsub("SWC_([^,]+,[^,]+)\\.csv", "\\1", basename(file))

  # Resetting the final 'DF' dataframe and variables
  DF[,2:45]<- NA
  rm(list = c("min_prev", "min_prev_val", "second_peak_doy", "second_peak_val", "first_peak_doy") [sapply(c("min_prev", "min_prev_val", "second_peak_doy", "second_peak_val", "first_peak_doy"), exists)])
  min_prev <- min_prev_val <- second_peak_doy <- second_peak_val <- first_peak_doy <- NA
  
  # seq_years <- 1960:1963
  # year_to_check <- 1971
  
  for (year_to_check in seq_years) {
    
    # Resetting the variables
    rm(list = c("minSOLDI", "min_B4_first","first_peak_doy", "first_peak_val","next_year_peak_doy", "next_peak_doy", "next_peak_val", "previous_peak_doy", "previous_peak_val" , "previous_year", "previous_year_peak_doy", "minDOY", "unimodal", "min_after_second_val", "min_after_second", "second_peak_doy", "second_peak_val")[sapply(c("minSOLDI", "min_B4_first","first_peak_doy", "first_peak_val","next_year_peak_doy", "next_peak_doy", "next_peak_val", "previous_peak_doy", "previous_peak_val" , "previous_year", "previous_year_peak_doy", "minDOY", "unimodal", "min_after_second_val","min_after_second", "second_peak_doy", "second_peak_val"), exists)])
    low_onset <- mod_onset <- extreme_onset <- low_offset <- mod_offset <- extreme_offset <- min_B4_first <- first_peak_doy <- first_peak_val <- next_year_peak_doy <- next_peak_doy <- next_peak_val <- previous_peak_doy <- previous_peak_val <- previous_year <- previous_year_peak_doy <- minDOY <- unimodal <- min_after_second_val <- min_after_second <- second_peak_doy <- second_peak_val <- minloess <- NA
    
    # Extracting the 50 first days of the next year
    next_year <- year_to_check + 1
    
    next_year_data <-res.model[res.model$year == next_year,]
    
    next_year_subset <- res.model %>%
      filter(year == next_year) %>%
      slice_head(n = 50)
    
    # Creating the dataframe containing the current year to proceed
    data_year <- res.model %>% 
      filter(year == year_to_check) %>% 
      mutate(DOY = row_number()) 
    data_year <- bind_rows(data_year, next_year_subset) %>%
      mutate(DOY = row_number())
    
    cumulativeSOLDI <- sum(data_year$SOLDI[data_year$year == year_to_check])
    DF$Cumulative.Year[DF$Year== year_to_check] <- cumulativeSOLDI
    
    # Smoothing curve parameters
    loess_model <- loess(SOLDI ~ DOY, data = data_year,
                         span = 0.16, degree = 1, family = "gaussian")
    
    fitted_values <- loess_model$fitted
    DOY<- c(1:length(fitted_values))
    df <- as.data.frame(cbind(DOY,fitted_values))
    
    # Peaks detection
    loess_peak <- peakfunction(df, threshold = 8, m = 20, max.peak = 3)
    
    peakA <- loess_peak[1,1]
    peakB <- loess_peak[2,1]
    peakC <- loess_peak[3,1]
    
    peakA_val <-  loess_peak[1,2]
    peakB_val <-  loess_peak[2,2]
    peakC_val <-  loess_peak[3,2]
    
    # Making sure the first peak isn't detected before previous year's last minimum value --> avoiding the first detected peak to be the same as the previous year's last peak.
    if (exists("min_prev") && !is.na(min_prev) && exists("peakA") && !is.na(peakA) && min_prev > peakA) {
      peakA <- peakB
      peakB <- peakC
      peakC <- NA
      
      peakA_val <- peakB_val
      peakB_val <- peakC_val
      peakC_val <- NA
    }  
    
    # Making sure a small early peak in beginning of the year isn't considered the first drought period.
    if (exists("peakA") && !is.na(peakA) && exists("peakB") && !is.na(peakB) && peakA < 100 && peakA_val/peakB_val < 0.8 & (min(df$fitted_values[df$DOY > peakA & df$DOY < peakB]))/peakA_val > 0.90 && (min(df$fitted_values[df$DOY < peakA]))/peakA_val > 0.90) {
      peakA <- peakB
      peakB <- peakC
      peakC <- NA
      
      peakA_val <- peakB_val
      peakB_val <- peakC_val
      peakC_val <- NA
    }
      
    # Making sure the first peak detected isn't the same as the previous year's last peak
    if (year_to_check != min(seq_years) && !is.na(peakA) && (!is.na(DF$First.Peak.DOY[DF$Year== year_to_check -1]) | !is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1])) ) {
      if (peakA < DF$First.Peak.DOY[DF$Year== year_to_check -1] - 365 + 20 | !is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]) && peakA < DF$Second.Peak.DOY[DF$Year== year_to_check -1] - 365 + 20){
        peakA <- peakB
        peakB <- peakC
        peakC <- NA
        
        peakA_val <- peakB_val
        peakB_val <- peakC_val
        peakC_val <- NA
      }
    }

    # Merging close peaks together --> making sure two peaks aren't considered two different droughts when they actually are the same drought
    if(!is.na(peakC)){
      if (min(peakA_val,peakB_val) / max(peakA_val,peakB_val) > 0.85 & min(peakB_val,peakC_val) / max(peakB_val,peakC_val) > 0.85 & (abs(min(df$fitted_values[df$DOY >= peakA & df$DOY <= peakB]))) / (abs(min(c(peakA_val, peakB_val)))) > 0.90 & (abs(min(df$fitted_values[df$DOY >= peakB & df$DOY <= peakC]))) / (abs(min(c(peakB_val, peakC_val)))) > 0.90) {
        min_B4_first <- which.min(data_year$SOLDI[data_year$year == year_to_check][0:peakA]) 
        min_B4_first_val <- data_year$SOLDI[data_year$year == year_to_check & data_year$DOY == min_B4_first]
        peakA <-  loess_peak[which.max(loess_peak[,2]), 1]
        fitted_values2 <- loess_model$fitted [peakC:length(data_year$DOY)]
        DOY2<- c(peakC:length(data_year$DOY))
        df2 <- as.data.frame(cbind(DOY2,fitted_values2))
        loess_peak2 <- peakfunction(df2, threshold = 8, m = 20, max.peak = 2)
        if (!is.na(loess_peak2[,2][1]) | !is.na(loess_peak2[,2][2])) {
          peakB <- loess_peak2[which.max(loess_peak2[,2]), 1]
          minloess <- min(df$fitted_values[][peakA:peakB])
          minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA & df$DOY < peakB ][1]
          if (abs(loess_peak2[loess_peak2[,1] == peakB,2][1] - minloess) < 0.9) {
            peakB <- NA
            unimodal <- TRUE
            minloess <- min(df$fitted_values[][peakA:length(df$fitted_values)])
            minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA][1]
          }
        } else {
          peakB <-  NA
          unimodal <- TRUE
          minloess <- min(df$fitted_values[][peakA:length(df$fitted_values)])
          minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA][1]
        }
        
        loess_peak[,1][1] <- peakA
        loess_peak[,1][2] <- peakB
        loess_peak[,1][3] <- NA
        
      } else if (peakB-peakA < 100 & min(peakA_val,peakB_val) / max(peakA_val,peakB_val) > 0.85 & (abs(min(df$fitted_values[df$DOY >= peakA & df$DOY <= peakB]))) / (abs(min(c(peakA_val, peakB_val)))) > 0.85 ) {
        minloess <- min(df$fitted_values[][peakB:peakC])
        minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakB][1]
        peakA <- data_year$DOY[data_year$DOY > (peakA - 20) & data_year$DOY < (peakB + 30) & data_year$SOLDI == max(data_year$SOLDI[data_year$DOY > (peakA - 20) & data_year$DOY < (peakB + 30)])][1]
        
      } else if (peakC-peakB < 100 & min(peakB_val,peakC_val) / max(peakB_val,peakC_val) > 0.85 & (abs(min(df$fitted_values[df$DOY >= peakB & df$DOY <= peakC]))) / (abs(min(c(peakB_val, peakC_val)))) > 0.85 ) {
        minloess <- min(df$fitted_values[][peakA:peakB])
        minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA][1]
        
      } else if (abs(peakC_val-peakA_val) < 3)  {
        minloess <- min(c(min(df$fitted_values[][peakA:peakB]),min(df$fitted_values[][peakB:peakC])))
        minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA][1]
        
      } else {
        loess_subset <- loess_peak[loess_peak[, 2] != min(loess_peak[, 2], na.rm = TRUE), ]
        
        peakA <- min(loess_subset[,1])
        peakA_val <- loess_subset[,2] [loess_subset[,1] == min(loess_subset[,1])]
        peakB <- max(loess_subset[,1])
        peakB_val <- loess_subset[,2] [loess_subset[,1] == max(loess_subset[,1])]
        minloess <- min(df$fitted_values[][peakA:peakB]) 
        minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA & df$DOY < peakB][1]
        min_B4_first <- min_B4_first_val <- NA
      }
    }
    
    
    if(!is.na(peakB) & is.na(peakC)){
      if (min(peakA_val,peakB_val) / max(peakA_val,peakB_val) > 0.85 & (abs(min(df$fitted_values[df$DOY >= peakA & df$DOY <= peakB]))) / (abs(min(c(peakA_val, peakB_val)))) > 0.85 ) {
        minloess <- NA
        unimodal <- TRUE
        
      } else minloess <- min(df$fitted_values[][peakA:peakB])
      minloess_DOY <- df$DOY[df$fitted_values == minloess & df$DOY > peakA & df$DOY < peakB][1]
      
    } 
    
    if(is.na(peakB)){
      minloess <- NA
      unimodal <- TRUE
    }  
    
    # When no peaks are detected (either ignoring the year when the drought index values are too low, or reprocessing with more sensitive LOESS parameters when the drought index values are sufficiently high)
    if(is.na(peakA) & is.na(peakB) & max(data_year$SOLDI[data_year$year == year_to_check]) > low_onset_threshold) {
      unimodal=TRUE
      
      loess_model <- loess(SOLDI ~ DOY, data = data_year, span = 0.09, degree = 1, family = "gaussian")
      
      fitted_values <- loess_model$fitted
      DOY <- c(1:length(fitted_values))
      df <- as.data.frame(cbind(DOY,fitted_values))
      loess_peak <- peakfunction(df, threshold = 2.25, m = 10, max.peak = 5)
      
      if (exists("min_prev") && !is.na(min_prev)){
        peakA <- loess_peak[loess_peak[,1] > min_prev, ][which.max(loess_peak[loess_peak[,1] > min_prev, 2]), 1]
      } else {peakA <- loess_peak[which.max(loess_peak[,2]), 1]}
      
      if (!length(peakA) > 0) {
        DF$Unimodal[DF$Year== year_to_check] <- NA
        DF$Cumulative.Year[DF$Year== year_to_check] <- sum(data_year$SOLDI[data_year$DOY > 0 & data_year$year == year_to_check])
        rm(list = c( "second_peak_doy", "second_peak_val", "first_peak_doy") [sapply(c("second_peak_doy", "second_peak_val", "first_peak_doy"), exists)])
        min_prev <- min_prev_val <- 0
        
        next
        }
      first_peak_val <- max(data_year$SOLDI[][data_year$DOY > peakA - 25 & data_year$DOY < peakA + 25]) # ajouter ici une condition pour que peak A > a DF year -1 peaks
      first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val]
      
      minSOLDI <- min(data_year$SOLDI[first_peak_doy:length(data_year$DOY)])
      minDOY <- tail(data_year$DOY[data_year$SOLDI == minSOLDI])[1]
      
      second_peak_val <- second_peak_doy <- NA
      
    } else if (is.na(peakA) && is.na(peakB) && max(data_year$SOLDI < low_onset_threshold)) {
      DF$Unimodal[DF$Year== year_to_check] <- NA
      DF$Cumulative.Year[DF$Year== year_to_check] <- sum(data_year$SOLDI[data_year$DOY > 0 & data_year$year == year_to_check])
      rm(list = c( "second_peak_doy", "second_peak_val", "first_peak_doy") [sapply(c( "second_peak_doy", "second_peak_val", "first_peak_doy"), exists)])
      min_prev <- min_prev_val <- 0
      
      next
    } else first_peak_doy <- first_peak_val <- NA
    
    
    # Ignoring the year when no peak reaches the minimum drought threshold (25)
    if(!is.na(peakA_val) && max(peakA_val, peakB_val, peakC_val, na.rm=TRUE) < 25){
      next
    } 
    
    # Making sure only one value of 'minloess_doy' exists
    if (exists("minloess_DOY") && length(minloess_DOY) >1) {
      minloess_DOY <- minloess_DOY[1]
    }
    
    # Determination of the second drought peak value and DOY
    if (!is.na(minloess)) {
      miin <- floor(minloess_DOY - 5)
      maax <- floor(min(max(loess_peak[,1], na.rm = TRUE) + 30, (nrow(data_year) - 1)))
      range_indices <- miin:maax
      second_peak_val <- max(data_year$SOLDI[][range_indices])
      second_peak_doy <- data_year$DOY[data_year$SOLDI == second_peak_val & data_year$DOY >= miin & data_year$DOY <= maax] [1]
    }
    
    # Determination of the first peak value and DOY in the case of bimodal and unimodal droughts
    if (exists("first_peak_doy") && is.na(first_peak_doy)) {
      if (!is.na(minloess)) {
        if (peakA > 30) {
          ifelse(year_to_check == min(seq_years), 
                 first_peak_val <- max(data_year$SOLDI[][(peakA - 30) : minloess_DOY]),
                 first_peak_val <- max(data_year$SOLDI[data_year$DOY > ifelse(is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]), 
                                                                           1, 
                                                                           ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check - 1]), 
                                                                                  DF$Second.Peak.DOY[DF$Year== year_to_check - 1] - length(res.model$year[res.model$year == year_to_check - 1]) + 15, 
                                                                                  DF$First.Peak.DOY[DF$Year== year_to_check - 1] - length(res.model$year[res.model$year == year_to_check - 1]) + 15)) & data_year$DOY > (peakA-30) & data_year$DOY < minloess_DOY])) #[(peakA-30):(peakA+30)]))
        } else { ifelse(year_to_check == min(seq_years),
                      first_peak_val <- max(data_year$SOLDI[1:minloess_DOY]),
                      first_peak_val <- max(data_year$SOLDI[data_year$DOY > ifelse(is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                1,
                                                                                ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check - 1]), 
                                                                                       DF$Second.Peak.DOY[DF$Year== year_to_check - 1] - length(res.model$year[res.model$year == year_to_check - 1]) + 15, 
                                                                                       DF$First.Peak.DOY[DF$Year== year_to_check - 1] - length(res.model$year[res.model$year == year_to_check - 1]) + 15)) & data_year$DOY > (peakA-30) & data_year$DOY < minloess_DOY & data_year$DOY > 1])) #][1:(peakA+30)]))
        }
        
        first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val  & data_year$DOY >= (peakA - 30) & data_year$DOY <= minloess_DOY] [1]
        
        minSOLDI <- min(data_year$SOLDI[first_peak_doy:second_peak_doy])
        minDOY <- data_year$DOY[data_year$SOLDI == minSOLDI] [1]
        unimodal <- FALSE
        
        if (minSOLDI == 0) {
          datemin <- data_year$DOY[data_year$SOLDI == 0 & data_year$DOY > first_peak_doy & data_year$DOY < second_peak_doy]
          min_datemin <- min(datemin)
          max_datemin <- max(datemin)
          median_datemin <- (min_datemin + max_datemin) / 2
          minDOY <- datemin[which.min(abs(datemin - median_datemin))]
        }
      } else {
        unimodal <- TRUE
        if (!is.na(peakB)) {
          first_peak_val <- max(data_year$SOLDI[][data_year$DOY > peakA - 40 & data_year$DOY < peakB + 40 & data_year$DOY > ifelse(year_to_check != min(seq_years),
                                                                                                                                ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                       DF$Second.Peak.DOY[DF$Year== year_to_check - 1] + 15 - 365,
                                                                                                                                       ifelse(!is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                              DF$First.Peak.DOY[DF$Year== year_to_check - 1] + 15 -365,
                                                                                                                                              0)),
                                                                                                                                0)])
         
          first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val & data_year$DOY > peakA - 40 & data_year$DOY < peakB + 40] [1]
          minloess_DOY <- NA
          minSOLDI <- min(data_year$SOLDI[data_year$DOY > first_peak_doy & data_year$DOY <= length(data_year$SOLDI)])
          minDOY <- data_year$DOY[data_year$SOLDI == minSOLDI] [1]
          
          min_B4_first <- which.min(data_year$SOLDI[data_year$DOY < first_peak_doy]) 
          min_B4_first_val <- data_year$SOLDI[data_year$year == year_to_check & data_year$DOY == min_B4_first]
        } else {
          first_peak_val <- max(data_year$SOLDI[data_year$SOLDI > c(100, data_year$SOLDI[-length(data_year$SOLDI)]) & data_year$DOY > peakA - 40 & data_year$DOY < peakA + 40 & data_year$DOY > ifelse(year_to_check != min(seq_years),
                                                                                                                              ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                     DF$Second.Peak.DOY[DF$Year== year_to_check - 1] + 15 - 365,
                                                                                                                                     ifelse(!is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                            DF$First.Peak.DOY[DF$Year== year_to_check - 1] + 15 -365,
                                                                                                                                            0)),
                                                                                                                              0)])
           c(NA, data_year$SOLDI[-length(data_year$SOLDI)])
          
          first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val  & data_year$DOY > peakA - 40 & data_year$DOY < peakA + 40] [1]
          minloess_DOY <- NA
          minSOLDI <- min(data_year$SOLDI[data_year$DOY > first_peak_doy & data_year$DOY <= length(data_year$SOLDI)])
          minDOY <- data_year$DOY[data_year$SOLDI == minSOLDI & data_year$DOY > first_peak_doy] [1]
          min_B4_first <- min_B4_first_val <- NA
          
          if (minSOLDI > first_peak_val) {
            first_peak_val <- max(data_year$SOLDI[data_year$DOY > peakA - 20 & data_year$DOY < peakA + 20 & data_year$DOY > ifelse(year_to_check != min(seq_years),
                                                                                                                                ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                       DF$Second.Peak.DOY[DF$Year== year_to_check - 1] + 15 - 365,
                                                                                                                                       ifelse(!is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]),
                                                                                                                                              DF$First.Peak.DOY[DF$Year== year_to_check - 1] + 15 -365,
                                                                                                                                              0)),
                                                                                                                                0)])
            first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val & data_year$DOY > peakA - 20 & data_year$DOY < peakA + 20] [1]
            minSOLDI <- min(data_year$SOLDI[data_year$DOY > first_peak_doy & data_year$DOY <= length(data_year$SOLDI)])
            minDOY <- data_year$DOY[data_year$SOLDI == minSOLDI & data_year$DOY > first_peak_doy] [1]

          }
        }
        
        if (minSOLDI == 0) {
          minDOY <- first_peak_doy + tail(which(data_year$SOLDI[first_peak_doy:length(data_year$SOLDI)] == 0), 1) - 1
        }
        second_peak_val <- NA
        second_peak_doy <- NA
        
      }

    }

    # In case the first drought peak is the same DOY as the minimum value
    if (first_peak_doy == minDOY) {
      first_peak_val <- max(data_year$SOLDI[data_year$DOY > peakA - 40 & data_year$DOY < peakA + (which.min(loess_model$fitted[peakA:(peakA + 40)])) - 1])
      first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val]
      minloess_DOY <- NA
      minSOLDI <- min(data_year$SOLDI[first_peak_doy:length(data_year$SOLDI)])
      minDOY <- tail(data_year$DOY[data_year$SOLDI == minSOLDI],1) [1]
      
    }
    
    # Determination of the min value and DOY before the first peak
    if (!exists("min_B4_first")) {
      min_B4_first <- which.min(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <  first_peak_doy]) 
      min_B4_first_val <- data_year$SOLDI[data_year$year == year_to_check & data_year$DOY == min_B4_first]
    } else if (exists("min_B4_first") & is.na(min_B4_first)) {
      min_B4_first <- which.min(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <  first_peak_doy]) 
      min_B4_first_val <- data_year$SOLDI[data_year$year == year_to_check & data_year$DOY == min_B4_first]
    }
    

    # Determination of the minimum value after the 2nd peak if there is any 2nd peak
    if (!is.na(second_peak_doy) && second_peak_doy >= 1 && second_peak_doy <= length(data_year$year)) {
      range_indices <- (second_peak_doy +1):length(data_year$year)
      
      if (length(range_indices) == 0) {
        min_after_second <- NA
        min_after_second_val <- NA
      } else {
        if (any(data_year$SOLDI[range_indices] == 0)) {
          min_after_second <- second_peak_doy  + (tail(which(data_year$SOLDI[range_indices] == 0), 1))
          min_after_second_val <- data_year$SOLDI[data_year$DOY == min_after_second]
        } else {
          min_after_second <- second_peak_doy + (which.min(data_year$SOLDI[range_indices]))
          min_after_second_val <- data_year$SOLDI[data_year$DOY == min_after_second]
        }
      }
    } else {
      min_after_second <- NA
      min_after_second_val <- NA
    }
    
    # In case the second drought peak is the same value as the minimum after the second drought
    if (!is.na(second_peak_val) && exists("min_after_second_val") && !is.na(min_after_second_val) & min_after_second_val >= second_peak_val) {
      miin <- floor(minloess_DOY - 5)
      maax <- floor(min(max(loess_peak[,1], na.rm = TRUE) + 30, (nrow(data_year) - 30)))
      range_indices <- miin:maax
      second_peak_val <- max(data_year$SOLDI[][range_indices])
      second_peak_doy <- data_year$DOY[data_year$SOLDI == second_peak_val  & data_year$DOY >= miin & data_year$DOY <= maax] [1]
      
      min_after_second_val <- min(data_year$SOLDI[data_year$DOY > second_peak_doy])
      min_after_second <- data_year$DOY[data_year$SOLDI == min_after_second_val]
      
    }
    
    # In case of a multi-year drought
      # Multi-year drought continuing after the current year

    if (!is.na(min_after_second) && min_after_second_val > low_offset_threshold & unimodal == FALSE) {
      multiyear_drought <-  TRUE
      
      data_year <- res.model %>%
        filter(between(year, year_to_check, year_to_check + 1 )) %>%
        mutate(DOY = row_number())
      
      # plot(data_year$SOLDI, type = 'l',
      #      main = paste("LOESS Tendency curve for year", year_to_check, "&", year_to_check+1),
      #      xlab = "DOY", ylab = "SOLDI")
      
      if (year_to_check == max(seq_years)) { 
        loess_model <- loess(SOLDI ~ DOY, data = data_year,
                             span = 0.16, degree = 1, family = "gaussian")
        
      } else {    
        loess_model <- loess(SOLDI ~ DOY, data = data_year,
                             span = 0.09, degree = 1, family = "gaussian")
      }
      
      fitted_values <- loess_model$fitted
      DOY<- c(1:length(fitted_values))
      df <- as.data.frame(cbind(DOY,fitted_values))
      loess_peak <- peakfunction(df, threshold = 2.25, m = 20, max.peak = 5)
      
      peakA <- ifelse(is.na(peakA),0,peakA )
      peakB <- ifelse(is.na(peakB),0,peakB )
      peakC <- ifelse(is.na(peakC),0,peakC )
      
      
      next_year_peak_doy <-  loess_peak[,1] [loess_peak[,1] > (min_after_second) & loess_peak[,1] > peakA & loess_peak[,1] > peakB & loess_peak[,1] > peakC ] [1]
      
      next_peak_val <- max(data_year$SOLDI[][data_year$DOY > next_year_peak_doy - 25 & data_year$DOY < next_year_peak_doy + 25 & data_year$DOY > min_after_second] )
      next_peak_doy <- data_year$DOY[data_year$SOLDI == next_peak_val & data_year$DOY > next_year_peak_doy - 25 & data_year$DOY < next_year_peak_doy + 25] [1]
      
      if (!is.na(next_year_peak_doy)) {
        min_after_second_val <- min(data_year$SOLDI[data_year$DOY >= min_after_second & data_year$DOY < next_peak_doy])
        min_after_second <- data_year$DOY[data_year$SOLDI == min_after_second_val & data_year$DOY > second_peak_doy][1]
      }
      
    } else if (!is.na(minDOY) && minSOLDI > low_offset_threshold & unimodal == TRUE) {
      data_year <- res.model %>%
        filter(between(year, year_to_check, year_to_check + 1 )) %>%
        mutate(DOY = row_number())
      
      if (year_to_check == max(seq_years)) { 
        loess_model <- loess(SOLDI ~ DOY, data = data_year, span = 0.16, degree = 1, family = "gaussian")
      } else {    
        loess_model <- loess(SOLDI ~ DOY, data = data_year, span = 0.09, degree = 1, family = "gaussian")
      }
      
      fitted_values <- loess_model$fitted
      DOY <- c(1:length(fitted_values))
      df <- as.data.frame(cbind(DOY,fitted_values))
      loess_peak <- peakfunction(df, threshold = 2.25, m = 20, max.peak = 5)
      
      peakA <- ifelse(is.na(peakA),0,peakA ) 
      peakB <- ifelse(is.na(peakB),0,peakB ) 
      peakC <- ifelse(is.na(peakC),0,peakC ) 
      
      next_year_peak_doy <-  loess_peak[,1] [loess_peak[,1] > minDOY & loess_peak[,1] > peakA & loess_peak[,1] > peakB & loess_peak[,1] > peakC & loess_peak[,1] > length(data_year$year[data_year$year == year_to_check])] [1]
      next_peak_val <- max(data_year$SOLDI[][data_year$DOY > minDOY & data_year$DOY > next_year_peak_doy - 25 & data_year$DOY < next_year_peak_doy + 25])
      next_peak_doy <- data_year$DOY[data_year$SOLDI == next_peak_val & data_year$DOY > next_year_peak_doy - 30 & data_year$DOY < next_year_peak_doy + 30] [1]
      
      if (!is.na(next_year_peak_doy)) {
        minSOLDI <- min(data_year$SOLDI[data_year$DOY > first_peak_doy & data_year$DOY < next_peak_doy])
        minDOY <- tail(data_year$DOY[data_year$SOLDI == minSOLDI & data_year$DOY > first_peak_doy & data_year$DOY < next_peak_doy],1)
      }
    }

    # 
    # if (is.numeric(min_B4_first)) {
    #   previous_peak_doy <- TRUE
    #   
    #   data_temp <- res.model %>%
    #     filter(between(year, year_to_check -1, year_to_check + 1) & (year != year_to_check + 1 | DOY <= 50)) %>%
    #     mutate(DOY = row_number())
    #   
    #   min_B4_first_val <- min(data_temp$SOLDI[data_temp$DOY > ifelse()])
    # }
    
    # print (min_B4_first)
    
    
    
     # Multiyear drought that started before the current year
    if (!is.na(min_B4_first_val) && min_B4_first_val > low_offset_threshold) {
      multiyear_drought <-  TRUE
      
      data_year <- res.model %>%
        filter(between(year, year_to_check -1, year_to_check + 1) & (year != year_to_check + 1 | DOY <= 50)) %>%
        mutate(DOY = row_number())
      
      loess_model <- loess(SOLDI ~ DOY, data = data_year,
                           span = 0.09, degree = 1, family = "gaussian")   
      
      # plot(data_year$SOLDI, type = 'l',
      #      main = paste("LOESS Tendency curve for year", year_to_check -1, "&", year_to_check),
      #      xlab = "DOY", ylab = "SOLDI")
      # lines(predict(loess_model, newdata = data_year), col = "red", lwd = 2)
      
      fitted_values <- loess_model$fitted
      DOY<- c(1:length(fitted_values))
      df <- as.data.frame(cbind(DOY,fitted_values))
      loess_peak <- peakfunction(df, threshold = 2.25, m = 20, max.peak = 5)
      
      first_peak_doy <- first_peak_doy + length(data_year$SOLDI [data_year$year == year_to_check-1])
      minDOY <- minDOY + length(data_year$SOLDI [data_year$year == year_to_check-1])

      if (exists("second_peak_doy") && !is.na(second_peak_doy)) {
        second_peak_doy <- second_peak_doy + length(data_year$SOLDI [data_year$year == year_to_check-1])
        min_after_second <- min_after_second + length(data_year$SOLDI [data_year$year == year_to_check-1])
      }
      previous_peak_doy <- ifelse(is.na(DF$First.Peak.DOY[DF$Year== year_to_check - 1]),
                                  1,
                                  max(DF$First.Peak.DOY[DF$Year== year_to_check-1],DF$Second.Peak.DOY[DF$Year== year_to_check-1], na.rm = TRUE))
        
      
      if (previous_peak_doy == first_peak_doy) {
        first_peak_val <- max(data_year$SOLDI[data_year$DOY > (which(loess_model$fitted[-length(loess_model$fitted)] < loess_model$fitted[-1])[ which(loess_model$fitted[-length(loess_model$fitted)] < loess_model$fitted[-1]) > 397][1]) & data_year$DOY < (which(loess_model$fitted[-length(loess_model$fitted)] < loess_model$fitted[-1])[ which(loess_model$fitted[-length(loess_model$fitted)] < loess_model$fitted[-1]) > 397] [1] + 60)])
        first_peak_doy <- data_year$DOY[data_year$SOLDI == first_peak_val]
      }
      
      min_B4_first_val <- min(data_year$SOLDI[data_year$DOY >= previous_peak_doy & data_year$DOY <= ifelse(previous_peak_doy > min_B4_first + length(data_year$year[data_year$year == year_to_check - 1]),
                                                                                                        first_peak_doy, 
                                                                                                        min_B4_first + length(data_year$year[data_year$year == year_to_check - 1])) ] , na.rm = TRUE)
      min_B4_first <-  tail(which(data_year$SOLDI ==  min_B4_first_val & data_year$DOY <= first_peak_doy),1) [1]
      
    } else { previous_peak_doy <- NA }
    
    # Making sure not to forget drought values between the previous year's last drought and the current year's first drought, if so, previous year's last drought offset phase is reprocessed
    if (year_to_check != min(seq_years) && exists("min_prev") && !is.na(min_prev)) {
      if (exists("previous_peak_doy") && !is.na(previous_peak_doy)){
        min_B4_first <- min_B4_first - length(data_year$year[data_year$year == year_to_check - 1])
      }
      if (min_B4_first < min_prev) {
        min_B4_first <- ifelse(min_prev > ifelse(exists("previous_peak_doy") & !is.na(previous_peak_doy), first_peak_doy - length(data_year$year[data_year$year == year_to_check - 1]), first_peak_doy), min_B4_first, min_prev)
        
      } else if (min_B4_first > min_prev){
        cumul <- (sum(res.model$SOLDI[res.model$year == year_to_check - 1 & res.model$DOY > (length(res.model$SOLDI[res.model$year == year_to_check]) + min_prev)]) +
                    sum(res.model$SOLDI[res.model$year == year_to_check & res.model$DOY < min_B4_first])) /(abs(min_prev) + abs(min_B4_first))
        if (cumul > 25 & min_prev_val > min_B4_first_val) {
          data_year_third <- res.model %>%
            filter(
              (year == year_to_check - 1 & DOY >= ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1])) |
                (year == year_to_check & DOY <= min_B4_first)
            ) %>%
            mutate(DOY = row_number())
          
          # plot(data_year_third$SOLDI, type = 'l',
          #      main = paste("LOESS Tendency curve for year", year_to_check -1, "&", year_to_check),
          #      xlab = "DOY", ylab = "SOLDI")
          
          peak <- which.max(data_year_third$SOLDI)
          peak_val <- data_year_third$SOLDI[data_year_third$DOY == peak]
          low_offset <- {tmp <- which(data_year_third$DOY > peak & data_year_third$SOLDI < low_offset_threshold); if (length(tmp) == 0) NA else tmp[1]}
          mod_offset <- {tmp <- which(data_year_third$DOY > peak & data_year_third$SOLDI < mod_offset_threshold); if (length(tmp) == 0) NA else tmp[1]}
          extreme_offset <- {tmp <- which(data_year_third$DOY > peak & data_year_third$SOLDI < high_offset_threshold); if (length(tmp) == 0) NA else tmp[1]}
          DS_low <- sum(data_year_third$SOLDI[data_year_third$SOLDI >= low_offset_threshold & data_year_third$DOY > length(data_year_third$DOY) - min_B4_first])
          DS_mod <- sum(data_year_third$SOLDI[data_year_third$SOLDI >= mod_offset_threshold & data_year_third$DOY > length(data_year_third$DOY) - min_B4_first])
          DS_extreme<- sum(data_year_third$SOLDI[data_year_third$SOLDI >= high_offset_threshold & data_year_third$DOY > length(data_year_third$DOY) - min_B4_first])
          
          if ( !is.na(DF$Unimodal[DF$Year== year_to_check - 1])){          
            if (DF$Unimodal[DF$Year== year_to_check - 1] == FALSE) {
              DF$Second.Low.D.S[DF$Year== year_to_check - 1] <- DF$Second.Low.D.S[DF$Year== year_to_check - 1] + DS_low
              DF$Second.Mod.D.S[DF$Year== year_to_check - 1] <- DF$Second.Mod.D.S[DF$Year== year_to_check - 1] + DS_mod
              DF$Second.Extreme.D.S[DF$Year== year_to_check - 1] <- DF$Second.Extreme.D.S[DF$Year== year_to_check - 1] + DS_extreme
              
              ifelse(DF$Second.Peak.Value[DF$Year== year_to_check - 1] >= low_offset_threshold, DF$Second.Low.D.Offset[DF$Year== year_to_check - 1] <- low_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Low.D.Offset[DF$Year== year_to_check - 1])
              ifelse(DF$Second.Peak.Value[DF$Year== year_to_check - 1] >= mod_offset_threshold, DF$Second.Mod.D.Offset[DF$Year== year_to_check - 1] <- mod_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Mod.D.Offset[DF$Year== year_to_check - 1])
              ifelse(DF$Second.Peak.Value[DF$Year== year_to_check - 1] >= high_offset_threshold, DF$Second.Extreme.D.Offset[DF$Year== year_to_check - 1] <- extreme_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Extreme.D.Offset[DF$Year== year_to_check - 1])
              
              DF$Second.Low.D.Duration[DF$Year== year_to_check - 1] <- DF$Second.Low.D.Offset[DF$Year== year_to_check - 1] - DF$Second.Low.D.Onset[DF$Year== year_to_check - 1]
              DF$Second.Mod.D.Duration[DF$Year== year_to_check - 1] <- DF$Second.Mod.D.Offset[DF$Year== year_to_check - 1] - DF$Second.Mod.D.Onset[DF$Year== year_to_check - 1]
              DF$Second.Extreme.D.Duration[DF$Year== year_to_check - 1] <- DF$Second.Extreme.D.Offset[DF$Year== year_to_check - 1] - DF$Second.Extreme.D.Onset[DF$Year== year_to_check - 1]
              
              DF$Second.Peak.DOY[DF$Year== year_to_check - 1] <- peak + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1
              DF$Second.Peak.Value[DF$Year== year_to_check - 1] <- peak_val
              
            } else {
              DF$First.Low.D.S[DF$Year== year_to_check - 1] <- DF$First.Low.D.S[DF$Year== year_to_check - 1] + DS_low
              DF$First.Mod.D.S[DF$Year== year_to_check - 1] <- DF$First.Mod.D.S[DF$Year== year_to_check - 1] + DS_mod
              DF$First.Extreme.D.S[DF$Year== year_to_check - 1] <- DF$First.Extreme.D.S[DF$Year== year_to_check - 1] + DS_extreme
              
              ifelse(DF$First.Peak.Value[DF$Year== year_to_check - 1] >= low_offset_threshold, DF$First.Low.D.Offset[DF$Year== year_to_check - 1] <- low_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Low.D.Offset[DF$Year== year_to_check - 1])
              ifelse(DF$First.Peak.Value[DF$Year== year_to_check - 1] >= mod_offset_threshold, DF$First.Mod.D.Offset[DF$Year== year_to_check - 1] <- mod_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Mod.D.Offset[DF$Year== year_to_check - 1])
              ifelse(DF$First.Peak.Value[DF$Year== year_to_check - 1] >= high_offset_threshold, DF$First.Extreme.D.Offset[DF$Year== year_to_check - 1] <- extreme_offset + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1, DF$First.Extreme.D.Offset[DF$Year== year_to_check - 1])
              
              DF$First.Low.D.Duration[DF$Year== year_to_check - 1] <- DF$First.Low.D.Offset[DF$Year== year_to_check - 1] - DF$First.Low.D.Onset[DF$Year== year_to_check - 1]
              DF$First.Mod.D.Duration[DF$Year== year_to_check - 1] <- DF$First.Mod.D.Offset[DF$Year== year_to_check - 1] - DF$First.Mod.D.Onset[DF$Year== year_to_check - 1]
              DF$First.Extreme.D.Duration[DF$Year== year_to_check - 1] <- DF$First.Extreme.D.Offset[DF$Year== year_to_check - 1] - DF$First.Extreme.D.Onset[DF$Year== year_to_check - 1]
              
              DF$First.Peak.DOY[DF$Year== year_to_check - 1] <- peak + ifelse(!is.na(DF$Second.Peak.DOY[DF$Year== year_to_check -1]), DF$Second.Peak.DOY[DF$Year== year_to_check -1], DF$First.Peak.DOY[DF$Year== year_to_check -1]) -1
              DF$First.Peak.Value[DF$Year== year_to_check - 1] <- peak_val
              
            }
          }
          
          low_offset <- mod_offset <- extreme_offset <- DS_low <- DS_mod <- DS_extreme <- NA
          
        } else if (cumul > 25 & min_prev_val < min_B4_first_val) {
          data_year <- res.model %>%
            filter(between(year, year_to_check -1, year_to_check + 1) & (year != year_to_check + 1 | DOY <= 50)) %>%
            mutate(DOY = row_number())
# 
#           plot(data_year$SOLDI, type = 'l',
#                main = paste("LOESS Tendency curve for year", year_to_check -1, "&", year_to_check),
#                xlab = "DOY", ylab = "SOLDI")
# 
#           first_peak_doy <- first_peak_doy + length(data_year$SOLDI [data_year$year == year_to_check-1])
#           minDOY <- minDOY + length(data_year$SOLDI [data_year$year == year_to_check-1])

          if (exists("previous_peak_doy") && is.na(previous_peak_doy)) {
            previous_peak_doy <- TRUE
            first_peak_doy <- first_peak_doy + length(data_year$SOLDI [data_year$year == year_to_check-1])
            minDOY <- minDOY + length(data_year$SOLDI [data_year$year == year_to_check-1])
            second_peak_doy <- second_peak_doy + length(data_year$SOLDI [data_year$year == year_to_check-1])
            min_after_second <- min_after_second + length(data_year$SOLDI [data_year$year == year_to_check-1])
          }
          min_B4_first <- min_prev
          min_B4_first_val <- min_prev_val
        }
      }
      if (exists("previous_peak_doy") && !is.na(previous_peak_doy)){
        min_B4_first <- min_B4_first + length(data_year$year[data_year$year == year_to_check-1])
      }
    }
    

    ### First drought onset phase computation ###
    # Dataframe containing values needed for the onset curve fitting (values before the minimum and after the peak are replaced by respectively the minimum and the peak values)
    data_year_first_drying <- data_year %>%
      mutate(SOLDI = ifelse(DOY > first_peak_doy, first_peak_val, SOLDI)) %>%
      mutate(SOLDI = ifelse(DOY <= min_B4_first, min_B4_first_val, SOLDI)) %>%
      mutate(SOLDI = SOLDI - min_B4_first_val) 
    
    # Extending the dataframe when the peak DOY is too late in the year for the curve to fit well
    if (data_year_first_drying$SOLDI[data_year_first_drying$DOY == max(data_year_first_drying$DOY) - 10] != data_year_first_drying$SOLDI[data_year_first_drying$DOY == max(data_year_first_drying$DOY)]){    
      new_rows <- data.frame(
        DOY = seq(length(data_year_first_drying$SOLDI)+1, length(data_year_first_drying$SOLDI) + 30), 
        SOLDI  = first_peak_val,
        date = seq(as.Date(max(data_year_first_drying$date),format="%Y-%m-%d")+1 ,as.Date(max(data_year_first_drying$date)+30,format="%Y-%m-%d"),by="day"),
        year = as.numeric(format(max(data_year_first_drying$date), "%Y"))
      )
    
      data_year_first_drying <- bind_rows(data_year_first_drying, new_rows)
    }
    
    # Elmore curve parameters
    M <- "Elmore"
    fit <- curvefit(data_year_first_drying$SOLDI, seq(1, length(data_year_first_drying$SOLDI), 1), methods = M)
    
    # l_pheno <- get_pheno(fit, M, IsPlot = FALSE)
    mn <- fit$model$Elmore$par[1]
    mx <- fit$model$Elmore$par[2]
    sos <- fit$model$Elmore$par[3]
    rsp <- fit$model$Elmore$par[4]
    eos <- fit$model$Elmore$par[5]
    rau <- fit$model$Elmore$par[6]
    m7 <- fit$model$Elmore$par[7]

    t <- seq(1, length(data_year_first_drying$SOLDI))
    
    elmore <- mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 1/(1 + exp(-rau * (t - eos))))
    
    drying_curve <- mn + (mx - m7 * t) * (1 / (1 + exp(-rsp * (t - sos))))
    drying_curve_der <- (-m7) * (1 / (1 + exp(-rsp * (t - sos)))) + (mx - m7 * t) * (rsp * exp(-rsp * (t - sos)) / (1 + exp(-rsp * (t - sos)))^2)
    
    mod_onset_from_der <- which(drying_curve_der == max(drying_curve_der, na.rm = TRUE))
    
    drying_rate <- drying_curve_der[mod_onset_from_der]
    
    if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
      intercept_drying_curve <- drying_curve[mod_onset_from_der] + drying_curve_der[mod_onset_from_der] * (365 - mod_onset_from_der) + min_B4_first_val
    } else {
      intercept_drying_curve <- drying_curve[mod_onset_from_der] - drying_curve_der[mod_onset_from_der] * mod_onset_from_der + min_B4_first_val
    }
    
    # ifelse(exists("previous_peak_doy") & !is.na(previous_peak_doy),
    #        intercept_drying_curve <- drying_curve[mod_onset_from_der] + drying_curve_der[mod_onset_from_der] * (365 - mod_onset_from_der) + min_B4_first_val,
    #        intercept_drying_curve <- drying_curve[mod_onset_from_der] - drying_curve_der[mod_onset_from_der] * mod_onset_from_der + min_B4_first_val)
    # 
    # 
    
    # plot(data_year_first_drying$DOY, data_year_first_drying$SOLDI,
    #      xlab = "DOY",
    #      ylab = "SOLDI",
    #      main = paste("SOLDI vs DOY (Year:", year_to_check, ")"),
    #      pch = 16,
    #      col = "black",
    #      type = "l")
    # 
    # lines(drying_curve + min_B4_first_val, type = 'l')
    # abline(v = mod_onset_from_der, col = "red", lty = 2)
  
    # Resetting the elmore and data_year_first_drying values
    elmore <- elmore + min_B4_first_val
    data_year_first_drying$SOLDI <-  data_year_first_drying$SOLDI + min_B4_first_val
    
    # Low drought onset DOY determination based on the fitting curve
    if (elmore[1] < low_onset_threshold & max(elmore) >= low_onset_threshold) {
      low_onset <- (which(elmore[min_B4_first:length(data_year$DOY)] > low_onset_threshold ) [1] ) + min_B4_first -1
      
      if (!is.na(low_onset)) {
        if (low_onset > first_peak_doy | first_peak_val < low_onset_threshold) {
          low_onset <- NA
        } 
      } 
    } else if (data_year_first_drying$SOLDI[1] > low_onset_threshold) {
      low_onset <- NA
    }
    
    # When the low onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
    if (is.na(low_onset) & data_year_first_drying$SOLDI[1] < low_onset_threshold & max(data_year_first_drying$SOLDI) > low_onset_threshold){
      low_onset <- which(data_year_first_drying$SOLDI >= low_onset_threshold) [1] 
    } else { low_onset <- low_onset
    }
    
    # Moderate drought onset DOY determination based on the fitting curve
    if (elmore[1] < mod_onset_threshold & max(elmore) >= mod_onset_threshold) {
      mod_onset <- (which(elmore[min_B4_first:length(data_year$DOY)] > mod_onset_threshold) [1]) + min_B4_first -1
      
      if (!is.na(mod_onset)) {
        if (mod_onset > first_peak_doy | first_peak_val < mod_onset_threshold) {
          mod_onset <- NA
        } 
      }
    } else if (data_year_first_drying$SOLDI[1] > mod_onset_threshold) {
      mod_onset <- NA
    }
    
    # When the moderate onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
    if (is.na(mod_onset) & data_year_first_drying$SOLDI[1] < mod_onset_threshold & max(data_year_first_drying$SOLDI) > mod_onset_threshold) {
      mod_onset <- which(data_year_first_drying$SOLDI >= mod_onset_threshold) [1]
    } else { mod_onset <- mod_onset
    }
    
    # Extreme drought onset DOY determination based on the fitting curve
    if (elmore[1] < high_onset_threshold & max(elmore) > high_onset_threshold) {
      extreme_onset <- (which(elmore[min_B4_first:length(data_year$DOY)] > high_onset_threshold) [1]) + min_B4_first -1
      
      if (!is.na(extreme_onset)) {
        if (extreme_onset > first_peak_doy | first_peak_val < high_onset_threshold) {
          extreme_onset <- NA
        } 
      } 
    }  else if (data_year_first_drying$SOLDI[1] > high_onset_threshold) {
      extreme_onset <- NA
    }
    
    # When the Extreme onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
    if (is.na(extreme_onset) & data_year_first_drying$SOLDI[1] < high_onset_threshold & max(data_year_first_drying$SOLDI) > high_onset_threshold){
      extreme_onset <- which(data_year_first_drying$SOLDI >= high_onset_threshold) [1]
    } else { extreme_onset <- extreme_onset
    }
    
    
    ### First drought offset phase computation ###
    # Extending the dataframe into the next year so the curve fits better when the minimum after the first peak occurs late
    if (unimodal == TRUE && minDOY > length(data_year$SOLDI)) {
      res_filtered <- res.model %>%
        filter(year == (year_to_check + 1), (DOY < minDOY - length(data_year$SOLDI) + 51), (DOY > 50) ) 
      
      data_year <- bind_rows(
        data_year,
        res_filtered
      ) %>%
        mutate(DOY = row_number())
    }
    
    # Dataframe containing values needed for the offset curve fitting (values after the minimum and before the peak are replaced by respectively the minimum and the peak values)
    data_year_first_wetting <- data_year %>%
      mutate(SOLDI = ifelse(DOY < first_peak_doy, first_peak_val, SOLDI)) %>%
      mutate(SOLDI = ifelse(DOY > minDOY, minSOLDI, SOLDI)) %>% 
      mutate(SOLDI = SOLDI - minSOLDI) 
    
    # Elmore curve parameters
    fit <- curvefit(data_year_first_wetting$SOLDI, seq(1, length(data_year_first_wetting$SOLDI), 1), methods = M)
    
    # l_pheno <- get_pheno(fit, M, IsPlot = FALSE)
    mn <- fit$model$Elmore$par[1]
    mx <- fit$model$Elmore$par[2]
    sos <- fit$model$Elmore$par[3]
    rsp <- fit$model$Elmore$par[4]
    eos <- fit$model$Elmore$par[5]
    rau <- fit$model$Elmore$par[6]
    m7 <- fit$model$Elmore$par[7]
    
    t <- seq(1, length(data_year_first_wetting$SOLDI))
    
    elmore <- mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 1/(1 + exp(-rau * (t - eos))))
    
    wetting_curve <- mn + (mx - m7 * t) * (1 - 1/(1 + exp(-rau * (t - eos))))
    wetting_curve_der <- (-m7) - (-m7 * (1 + exp(-rau * (t -  eos))) - (mx - m7 * t) * (-rau) * exp(-rau * (t - eos)))/(1 + exp(-rau * (t - eos)))^2
    
    mod_offset_from_der <- which(wetting_curve_der[] == min(wetting_curve_der[], na.rm = TRUE))
    
    wetting_rate <- wetting_curve_der[mod_offset_from_der]
    
    if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
      intercept_wetting_rate <- wetting_curve[mod_offset_from_der] + wetting_curve_der[mod_offset_from_der] * (365 - mod_offset_from_der) +  minSOLDI
    } else {
      intercept_wetting_rate <- wetting_curve[mod_offset_from_der] - wetting_curve_der[mod_offset_from_der] * mod_offset_from_der + minSOLDI
    }
    
    
    # ifelse(exists("previous_peak_doy") & !is.na(previous_peak_doy),
    #        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] + wetting_curve_der[mod_offset_from_der] * (365 - mod_offset_from_der) + minSOLDI,
    #        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] - wetting_curve_der[mod_offset_from_der] * mod_offset_from_der + minSOLDI)
    # 
    
    # plot(data_year_first_wetting$DOY, data_year_first_wetting$SOLDI,
    #      xlab = "DOY",
    #      ylab = "SOLDI",
    #      main = paste("SOLDI vs DOY (Year:", year_to_check, ")"),
    #      pch = 16,
    #      col = "blue",
    #      type = "l")
    # 
    # points(wetting_curve + minSOLDI, type = 'l')
    # abline(v = mod_offset_from_der, col = "red", lty = 2)
    
    # Resetting the elmore and data_year_first_wetting values
    elmore <- elmore + minSOLDI
    data_year_first_wetting$SOLDI <-  data_year_first_wetting$SOLDI + minSOLDI
    
    # Low drought offset DOY determination based on the fitting curve
    if (elmore[length(data_year$DOY)] > low_offset_threshold) {
      low_offset <- NA
    } else if (elmore[length(data_year$DOY)] <= low_offset_threshold & max(elmore) >= low_onset_threshold) {
      low_offsets_in_fitperiod <- which(elmore[1:minDOY] >= low_offset_threshold)
      low_offset <- tail(low_offsets_in_fitperiod, 1)
      if(low_offset < first_peak_doy | first_peak_val < low_offset_threshold){
        low_offset <- NA
      }
    }
    
    # When the low offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
    if (is.na(low_offset) & tail(data_year_first_wetting$SOLDI,1) < low_offset_threshold & max(data_year_first_wetting$SOLDI) > low_offset_threshold){
      low_offsets_in_fitperiod <- which(data_year_first_wetting$SOLDI >= low_offset_threshold)
      if (length(low_offsets_in_fitperiod) > 0) {
        
        low_offset <- tail(low_offsets_in_fitperiod,1)
      } else { low_offset <- NA
      }
    } else { low_offset <- low_offset
    }
    
    # Moderate drought onset DOY determination based on the fitting curve
    if (elmore[length(data_year$DOY)] > mod_offset_threshold) {
      mod_offset <- NA
    } else if (elmore[length(data_year$DOY)] <= mod_offset_threshold & max(elmore) >= mod_onset_threshold) {
      mod_offsets_in_fitperiod <- which(elmore[1:minDOY] >= mod_offset_threshold)
      mod_offset <- tail(mod_offsets_in_fitperiod,1)
      if(mod_offset < first_peak_doy | first_peak_val < mod_offset_threshold){
        mod_offset <- NA
      }
    }

    # When the moderate offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
    if (is.na(mod_offset) & tail(data_year_first_wetting$SOLDI,1) < mod_offset_threshold & max(data_year_first_wetting$SOLDI) > mod_offset_threshold){
      mod_offsets_in_fitperiod <- which(data_year_first_wetting$SOLDI >= mod_offset_threshold)
      if (length(mod_offsets_in_fitperiod) > 0) {
        
        mod_offset <- tail(mod_offsets_in_fitperiod,1)
      } else { mod_offset <- NA
      }
    } else { mod_offset <- mod_offset
    }
    
    # Extreme drought onset DOY determination based on the fitting curve
    if (elmore[length(data_year$DOY)] > high_offset_threshold) {
      extreme_offset <- NA
    } else if (elmore[length(data_year$DOY)] <= high_offset_threshold & max(elmore) >= high_onset_threshold) {
      extreme_offsets_in_fitperiod <- which(elmore[1:minDOY] >= high_offset_threshold)
      extreme_offset <- tail(extreme_offsets_in_fitperiod,1)
      if(extreme_offset < first_peak_doy | first_peak_val < high_offset_threshold){
        extreme_offset <- NA
      }
    }
    
    # When the extreme offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
    if (is.na(extreme_offset) & tail(data_year_first_wetting$SOLDI,1) < high_offset_threshold & max(data_year_first_wetting$SOLDI) > high_offset_threshold){
      extreme_offsets_in_fitperiod <- which(data_year_first_wetting$SOLDI >= high_offset_threshold)
      if (length(extreme_offsets_in_fitperiod) > 0) {
        
        extreme_offset <- tail(extreme_offsets_in_fitperiod,1)
      } else { extreme_offset <- NA
      }
    } else { extreme_offset <- extreme_offset
    }
    
    # Inserting first drought onsets and offsets into the final "DF" dataframe
    first_low_onset <- low_onset
    first_low_offset <- low_offset
    first_mod_onset <- mod_onset
    first_mod_offset <- mod_offset
    first_extreme_onset <- extreme_onset
    first_extreme_offset <- extreme_offset
    
    DF$First.Drying.Rate[DF$Year== year_to_check] <- drying_rate
    DF$First.Wetting.Rate[DF$Year== year_to_check] <- wetting_rate
    DF$First.Intercept.drying[DF$Year== year_to_check] <- intercept_drying_curve
    DF$First.Intercept.wetting[DF$Year== year_to_check] <- intercept_wetting_rate
    
    low_onset <- mod_onset <- extreme_onset <- low_offset <- mod_offset <- extreme_offset <-  drying_rate <- wetting_rate <- intercept_drying_curve <- intercept_wetting_rate <-NA
    
    # Only processing when the year is bimodal
    if (unimodal == FALSE) {
      ### Second drought onset phase computation ###
      # Dataframe containing values needed for the onset curve fitting (values before the minimum and after the peak are replaced by respectively the minimum and the peak values)
      data_year_second_drying <- data_year %>%
        mutate(SOLDI = ifelse(DOY < minDOY, minSOLDI, SOLDI)) %>%
        mutate(SOLDI = ifelse(DOY > second_peak_doy, second_peak_val, SOLDI)) %>%
        mutate(SOLDI = SOLDI - minSOLDI)

      # Elmore curve parameters
      fit <- curvefit(data_year_second_drying$SOLDI, seq(1, length(data_year_second_drying$SOLDI), 1), methods = M)

      #l_pheno <- get_pheno(fit, M, IsPlot = FALSE)
      mn <- fit$model$Elmore$par[1]
      mx <- fit$model$Elmore$par[2]
      sos <- fit$model$Elmore$par[3]
      rsp <- fit$model$Elmore$par[4]
      eos <- fit$model$Elmore$par[5]
      rau <- fit$model$Elmore$par[6]
      m7 <- fit$model$Elmore$par[7]

      t <- seq(1, length(data_year_second_drying$SOLDI))
      elmore <- mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 1/(1 + exp(-rau * (t - eos))))
      drying_curve <- mn + (mx - m7 * t) * (1 / (1 + exp(-rsp * (t - sos))))
      drying_curve_der <- (-m7) * (1 / (1 + exp(-rsp * (t - sos)))) + (mx - m7 * t) * (rsp * exp(-rsp * (t - sos)) / (1 + exp(-rsp * (t - sos)))^2)
      
      if (!any(is.finite(drying_curve_der))) {
        t <- seq(1, length(t) + ((sos - length(t))/2))
        drying_curve <- mn + (mx - m7 * t) * (1 / (1 + exp(-rsp * (t - sos))))
        drying_curve_der <- (-m7) * (1 / (1 + exp(-rsp * (t - sos)))) + (mx - m7 * t) * (rsp * exp(-rsp * (t - sos)) / (1 + exp(-rsp * (t - sos)))^2)
      }
      
      mod_onset_from_der <- which(drying_curve_der == max(drying_curve_der, na.rm = TRUE))
      drying_rate <- drying_curve_der[mod_onset_from_der]

      if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
        intercept_drying_curve <- drying_curve[mod_onset_from_der] +
          drying_curve_der[mod_onset_from_der] * (365 - mod_onset_from_der) + minSOLDI
      } else {
        intercept_drying_curve <- drying_curve[mod_onset_from_der] -
          drying_curve_der[mod_onset_from_der] * mod_onset_from_der + minSOLDI
      }
      
      # ifelse(exists("previous_peak_doy") & !is.na(previous_peak_doy),
      #        intercept_drying_curve <- drying_curve[mod_onset_from_der] + drying_curve_der[mod_onset_from_der] * (365 - mod_onset_from_der) + minSOLDI,
      #        intercept_drying_curve <- drying_curve[mod_onset_from_der] - drying_curve_der[mod_onset_from_der] * mod_onset_from_der + minSOLDI)

      # plot(data_year_second_drying$DOY, data_year_second_drying$SOLDI,
      #      xlab = "DOY",
      #      ylab = "SOLDI",
      #      main = paste("SOLDI vs DOY (Year:", year_to_check, ")"),
      #      pch = 16,
      #      col = "blue",
      #      type = "l")
      # 
      # points(drying_curve , type = 'l')
      # abline(v = mod_onset_from_der, col = "red", lty = 2)

      # Resetting the elmore and data_year_second_drying values
      elmore <- elmore + minSOLDI
      data_year_second_drying$SOLDI <-  data_year_second_drying$SOLDI + minSOLDI
      
      # Low drought onset DOY determination based on the fitting curve
      if (elmore[1] < low_onset_threshold & max(elmore) > low_onset_threshold) {
        low_onset <- (which(elmore[minDOY:length(data_year$DOY)] > low_onset_threshold) [1]) + minDOY -1

        if (!is.na(low_onset)) {
          if (low_onset > second_peak_doy | second_peak_val < low_onset_threshold) {
            low_onset <- NA
          }
        }
      }

      # When the low onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
      if (is.na(low_onset) & data_year_second_drying$SOLDI[1] < low_onset_threshold & max(data_year_second_drying$SOLDI) > low_onset_threshold){
        low_onset <- which(data_year_second_drying$SOLDI >= low_onset_threshold) [1]
      }


      # Moderate drought onset DOY determination based on the fitting curve
      if (elmore[1] < mod_onset_threshold & max(elmore)> mod_onset_threshold) {
        mod_onset <- (which(elmore[minDOY:length(data_year$DOY)] >= mod_onset_threshold) [1] ) + minDOY - 1

        if (!is.na(mod_onset)) {
          if (mod_onset > second_peak_doy | second_peak_val < mod_onset_threshold) {
            mod_onset <- NA
          }
        }
      }

      # When the moderate onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
      if (is.na(mod_onset) & data_year_second_drying$SOLDI[1] < mod_onset_threshold & max(data_year_second_drying$SOLDI) > mod_onset_threshold){
        mod_onset <- which(data_year_second_drying$SOLDI >= mod_onset_threshold)[1]
      }


      # Extreme drought onset DOY determination based on the fitting curve
      if (elmore[1] < high_onset_threshold & max(elmore)> high_onset_threshold) {
        extreme_onset <- (which(elmore[minDOY:length(data_year$DOY)] >= high_onset_threshold) [1] ) + minDOY -1

        if (!is.na(extreme_onset)) {
          if (extreme_onset > second_peak_doy | second_peak_val < high_onset_threshold) {
            extreme_onset <- NA
          }
        }
      }

      # When the extreme onset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY after the peak), determination of the onset based on the raw KBDI values
      if (is.na(extreme_onset) & data_year_second_drying$SOLDI[1] < high_onset_threshold & max(data_year_second_drying$SOLDI) > high_onset_threshold){
        extreme_onset <- which(data_year_second_drying$SOLDI >= high_onset_threshold) [1]
      }

      
      
      ### Second drought offset phase computation ###
      # Extending the dataframe into the next year so the curve fits better when the minimum after the first peak occurs late
      if (min_after_second >= length(data_year$SOLDI)) {
        res_filtered <- res.model %>%
          filter(year == (year_to_check + 1), (DOY < min_after_second - length(data_year$SOLDI) + 51), (DOY > 50) ) 
        
        data_year <- bind_rows(
          data_year,
          res_filtered
        ) %>%
          mutate(DOY = row_number())
      }

      # Dataframe containing values needed for the offset curve fitting (values after the minimum and before the peak are replaced by respectively the minimum and the peak values)
      data_year_second_wetting <- data_year %>%
        mutate(SOLDI = ifelse(DOY < second_peak_doy, second_peak_val, SOLDI)) %>%
        mutate(SOLDI = ifelse(DOY >= min_after_second, min_after_second_val, SOLDI)) %>%
        mutate(SOLDI = SOLDI - min_after_second_val)
      
      # if (year_to_check == max(seq_years)) {
      #   fit <- curvefit(data_year_second_wetting$SOLDI, seq(1, length(data_year_second_wetting$SOLDI) -1, 1), methods = M)
      #   
      # } else {      
      #   fit <- curvefit(data_year_second_wetting$SOLDI, seq(1, length(data_year_second_wetting$SOLDI), 1), methods = M)
      #   }

      # Elmore curve parameters
      fit <- curvefit(data_year_second_wetting$SOLDI, seq(1, length(data_year_second_wetting$SOLDI), 1), methods = M)
      
      #l_pheno <- get_pheno(fit, M, IsPlot = FALSE) 
      mn <- fit$model$Elmore$par[1]
      mx <- fit$model$Elmore$par[2]
      sos <- fit$model$Elmore$par[3]
      rsp <- fit$model$Elmore$par[4]
      eos <- fit$model$Elmore$par[5]
      rau <- fit$model$Elmore$par[6]
      m7 <- fit$model$Elmore$par[7]

      #
      if (eos < 0 | rau < 0) {
        fit <- curvefit(data_year_second_wetting$SOLDI, seq(1, length(data_year_second_wetting$SOLDI)-1, 1), methods = M)
        mn <- fit$model$Elmore$par[1]
        mx <- fit$model$Elmore$par[2]
        sos <- fit$model$Elmore$par[3]
        rsp <- fit$model$Elmore$par[4]
        eos <- fit$model$Elmore$par[5]
        rau <- fit$model$Elmore$par[6]
        m7 <- fit$model$Elmore$par[7] 
      }
      
      t <- seq(1, length(data_year_second_wetting$SOLDI))

      elmore <- mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 1/(1 + exp(-rau * (t - eos))))

      wetting_curve <- mn + (mx - m7 * t) * (1 - 1/(1 + exp(-rau * (t - eos))))
      wetting_curve_der <- (-m7) - (-m7 * (1 + exp(-rau * (t -  eos))) - (mx - m7 * t) * (-rau) * exp(-rau * (t - eos)))/(1 + exp(-rau * (t - eos)))^2

      mod_offset_from_der <- which(wetting_curve_der[] == min(wetting_curve_der[], na.rm = TRUE))

      wetting_rate <- wetting_curve_der[mod_offset_from_der]
      
      # Case when the Elmore curve doesn"t fit properly
      if (!any(is.finite(wetting_curve_der))){
        
      loess_model <- loess(SOLDI ~ DOY, data = data_year_second_wetting, span = 0.3, degree = 1, family = "gaussian")
      
      fit <- curvefit(predict(loess_model), seq(1, length(data_year_second_wetting$SOLDI), 1), methods = M)
      
      #l_pheno <- get_pheno(fit, M, IsPlot = FALSE) 
      mn <- fit$model$Elmore$par[1]
      mx <- fit$model$Elmore$par[2]
      sos <- fit$model$Elmore$par[3]
      rsp <- fit$model$Elmore$par[4]
      eos <- fit$model$Elmore$par[5]
      rau <- fit$model$Elmore$par[6]
      m7 <- fit$model$Elmore$par[7]
      
      t <- seq(1, length(data_year_second_wetting$SOLDI))
      
      elmore <- mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 1/(1 + exp(-rau * (t - eos))))
      
      wetting_curve <- mn + (mx - m7 * t) * (1 - 1/(1 + exp(-rau * (t - eos))))
      wetting_curve_der <- (-m7) - (-m7 * (1 + exp(-rau * (t -  eos))) - (mx - m7 * t) * (-rau) * exp(-rau * (t - eos)))/(1 + exp(-rau * (t - eos)))^2
      
      mod_offset_from_der <- which(wetting_curve_der[] == min(wetting_curve_der[], na.rm = TRUE))
      
      wetting_rate <- wetting_curve_der[mod_offset_from_der]
      
      lines(predict(loess_model, newdata = data_year), col = "red", lwd = 2)
      }
      
      # !!!
      if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] +
          wetting_curve_der[mod_offset_from_der] * (365 - mod_offset_from_der) + min_after_second_val
      } else {
        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] -
          wetting_curve_der[mod_offset_from_der] * mod_offset_from_der + min_after_second_val
      }
      
      # ifelse(exists("previous_peak_doy" & !is.na(previous_peak_doy)),
      #        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] + wetting_curve_der[mod_offset_from_der] * (365 - mod_offset_from_der) + min_after_second_val,
      #        intercept_wetting_rate <- wetting_curve[mod_offset_from_der] - wetting_curve_der[mod_offset_from_der] * mod_offset_from_der + min_after_second_val)

      # plot(data_year_second_wetting$DOY, data_year_second_wetting$SOLDI,
      #      xlab = "DOY",
      #      ylab = "SOLDI",
      #      main = paste("SOLDI vs DOY (Year:", year_to_check, ")"),
      #      pch = 16,
      #      col = "blue",
      #      type = "l")
      # 
      # points(wetting_curve, type = 'l')
      # abline(v = mod_offset_from_der, col = "red", lty = 2)

      # Resetting the elmore and data_year_second_wetting values
      elmore <- elmore + min_after_second_val
      data_year_second_wetting$SOLDI <-  data_year_second_wetting$SOLDI + min_after_second_val

      # Low drought offset DOY determination based on the fitting curve
      if (elmore[length(data_year$DOY)] > low_offset_threshold) {
        low_offset <- NA
      } else if (elmore[length(data_year$DOY)] <= low_offset_threshold & max(data_year_second_wetting$SOLDI) >= low_offset_threshold) {
        low_offsets_in_fitperiod <- which(elmore[1:min_after_second] >= low_offset_threshold)
        low_offset <- if (length(low_offsets_in_fitperiod) > 0) tail(low_offsets_in_fitperiod, 1) else NA
        if(low_offset < second_peak_doy | second_peak_val < low_offset_threshold){
          low_offset <- NA
        }
      }

      # When the low offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
      if (is.na(low_offset) & tail(data_year_second_wetting$SOLDI,1) < low_offset_threshold & max(data_year_second_wetting$SOLDI) > low_offset_threshold){
        low_offsets_in_fitperiod <- which(data_year_second_wetting$SOLDI[1:length(data_year$DOY)] >= low_offset_threshold)
        if (length(low_offsets_in_fitperiod) > 0) {

          low_offset <- tail(low_offsets_in_fitperiod,1)
        } else { low_offset <- NA
        }
      } else { low_offset <- low_offset
      }


      # Moderate drought offset DOY determination based on the fitting curve
      if (elmore[length(data_year$DOY)] > mod_offset_threshold) {
        mod_offset <- NA
      } else if (elmore[length(data_year$DOY)] <= mod_offset_threshold & max(data_year_second_wetting$SOLDI) >= mod_offset_threshold) {
        mod_offsets_in_fitperiod <- which(elmore[1:min_after_second] >= mod_offset_threshold)
        mod_offset <- if (length(mod_offsets_in_fitperiod) > 0) tail(mod_offsets_in_fitperiod, 1) else NA
        if(mod_offset < second_peak_doy | second_peak_val < mod_offset_threshold){
          mod_offset <- NA
        }
      }

      # When the moderate offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
      if (is.na(mod_offset) & tail(data_year_second_wetting$SOLDI,1) < mod_offset_threshold & max(data_year_second_wetting$SOLDI) > mod_offset_threshold){
        mod_offsets_in_fitperiod <- which(data_year_second_wetting$SOLDI[1:length(data_year$DOY)] >= mod_offset_threshold)
        if (length(mod_offsets_in_fitperiod) > 0) {

          mod_offset <- tail(mod_offsets_in_fitperiod,1)
        } else { mod_offset <- NA
        }
      } else { mod_offset <- mod_offset
      }

      # Extreme drought offset DOY determination based on the fitting curve
      if (elmore[length(data_year$DOY)] > high_offset_threshold) {
        extreme_offset <- NA
      } else if (elmore[length(data_year$DOY)] <= high_offset_threshold & max(data_year_second_wetting$SOLDI) >= high_offset_threshold) {
        extreme_offsets_in_fitperiod <- which(elmore[1:min_after_second] >= high_offset_threshold)
        extreme_offset <- if (length(extreme_offsets_in_fitperiod) > 0) tail(extreme_offsets_in_fitperiod, 1) else NA
        if (!is.na(extreme_offset)){
          if (extreme_offset < second_peak_doy | second_peak_val < high_offset_threshold){
            extreme_offset <- NA
          }
        }
      }

      # When the extreme offset DOY hasn't been determined previously (either because the value is higher than the peak or the DOY before the peak), determination of the offset based on the raw KBDI values
      if (is.na(extreme_offset) & tail(data_year_second_wetting$SOLDI,1) < high_offset_threshold & max(data_year_second_wetting$SOLDI) > high_offset_threshold) {
        extreme_offsets_in_fitperiod <- which(data_year_second_wetting$SOLDI[1:length(data_year$DOY)] >= high_offset_threshold)
        if (length(extreme_offsets_in_fitperiod) > 0) {

          extreme_offset <- tail(extreme_offsets_in_fitperiod,1)
        } else { extreme_offset <- NA
        }
      } else {extreme_offset <- extreme_offset}
      
    
      # Wetting phase variables
      if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
        DF$Peak.Monsoon.DOY[DF$Year == year_to_check] <- minDOY - length(data_year$DOY[data_year$year == year_to_check - 1])
        DF$Wetting.Phase.Duration[DF$Year == year_to_check] <- (minDOY - length(data_year$DOY[data_year$year == year_to_check - 1])) - (first_peak_doy - length(data_year$DOY[data_year$year == year_to_check - 1]))
      } else {
        DF$Peak.Monsoon.DOY[DF$Year == year_to_check] <- minDOY
        DF$Wetting.Phase.Duration[DF$Year == year_to_check] <- minDOY - first_peak_doy
      }
      
      DF$Peak.Monsoon[DF$Year == year_to_check] <- first_peak_val - minSOLDI

  
      second_low_onset <- low_onset
      second_low_offset <- low_offset
      second_mod_onset <- mod_onset
      second_mod_offset <- mod_offset
      second_extreme_onset <- extreme_onset
      second_extreme_offset <- extreme_offset
    
    } else {
      second_low_onset <- low_onset  
      second_low_offset <- low_offset
      second_mod_onset <- mod_onset
      second_mod_offset <- mod_offset
      second_extreme_onset <- extreme_onset
      second_extreme_offset <- extreme_offset
      
      
    }

    # First low drought duration determination
    if (!is.na(first_low_onset)){
      if (!is.na(first_low_offset)) {
        first_low_duration <- first_low_offset - first_low_onset + 1
      } else if (exists("second_low_offset") && !is.na(second_low_offset)) { 
        first_low_duration <- second_low_offset - first_low_onset + 1
      } else first_low_duration <- NA
    } else first_low_duration <- NA
    
    # First moderate drought duration determination
    if (!is.na(first_mod_onset)){
      if (!is.na(first_mod_offset)) {
        first_mod_duration <- first_mod_offset - first_mod_onset + 1
      } else if (exists("second_mod_offset") && !is.na(second_mod_offset)) { 
        first_mod_duration <- second_mod_offset - first_mod_onset + 1
      } else first_mod_duration <- NA
    } else first_mod_duration <- NA

    # First extreme drought duration determination
    if (!is.na(first_extreme_onset)){
      if (!is.na(first_extreme_offset)) {
        first_extreme_duration <- first_extreme_offset - first_extreme_onset + 1
      } else if (exists("second_extreme_offset") && !is.na(second_extreme_offset)) { 
        first_extreme_duration <- second_extreme_offset - first_extreme_onset + 1
      } else first_extreme_duration <- NA
    } else first_extreme_duration <- NA
    
    # Second low drought duration determination
    if (exists("second_low_offset") && !is.na(second_low_offset)){
      if (!is.na(second_low_onset)) {
        second_low_duration <- second_low_offset - second_low_onset + 1
      } else second_low_duration <- second_low_offset - first_low_onset + 1
    } else second_low_duration <- NA 
    
    # Second moderate drought duration determination
    if (exists("second_mod_offset") &&!is.na(second_mod_offset)){
      if (!is.na(second_mod_onset)) {
        second_mod_duration <- second_mod_offset - second_mod_onset + 1
      } else second_mod_duration <- second_mod_offset - first_mod_onset + 1
    } else second_mod_duration <- NA
    
    # Second extreme drought duration determination
    if (exists("second_extreme_offset") && !is.na(second_extreme_offset)){
      if (!is.na(second_extreme_onset)) {
        second_extreme_duration <- second_extreme_offset - second_extreme_onset + 1
      } else second_extreme_duration <- second_extreme_offset - first_extreme_onset + 1
    } else second_extreme_duration <- NA
    
    # First low drought severity determination
    if (!is.na(first_low_onset)){
      if (!is.na(first_low_offset)) {
        DS_first_low <- sum(data_year$SOLDI[first_low_onset:first_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (exists("second_low_offset") && !is.na(second_low_offset)) { 
        DS_first_low <- sum(data_year$SOLDI[first_low_onset:second_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_first_low <- sum(data_year$SOLDI[first_low_onset:length(data_year$DOY[data_year$year == year_to_check])]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_low_offset)) {
      DS_first_low <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= first_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_peak_val) && first_peak_val > low_onset_threshold) {
      DS_first_low <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_first_low <- NA
    
    # First moderate drought severity determination
    if (!is.na(first_mod_onset)){
      if (!is.na(first_mod_offset)) {
        DS_first_mod <- sum(data_year$SOLDI[first_mod_onset:first_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (exists("second_mod_offset") && !is.na(second_mod_offset)) { 
        DS_first_mod <- sum(data_year$SOLDI[first_mod_onset:second_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_first_mod <- sum(data_year$SOLDI[first_mod_onset:length(data_year$DOY[data_year$year == year_to_check])]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_mod_offset)) {
      DS_first_mod <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= first_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_peak_val) && first_peak_val > mod_onset_threshold) {
      DS_first_mod <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_first_mod <- NA
    
    # First extreme drought severity determination
    if (!is.na(first_extreme_onset)){
      if (!is.na(first_extreme_offset)) {
        DS_first_extreme <- sum(data_year$SOLDI[first_extreme_onset:first_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (exists("second_extreme_offset") && !is.na(second_extreme_offset)) { 
        DS_first_extreme <- sum(data_year$SOLDI[first_extreme_onset:second_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_first_extreme <- sum(data_year$SOLDI[first_extreme_onset:length(data_year$DOY[data_year$year == year_to_check])]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_extreme_offset)) {
      DS_first_extreme <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= first_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(first_peak_val) && first_peak_val > high_onset_threshold) {
      DS_first_extreme <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_first_extreme <- NA
    
    # Second low drought severity determination
    if (exists("second_low_offset") && !is.na(second_low_offset)){
      if (!is.na(second_low_onset)) {
        DS_second_low <- sum(data_year$SOLDI[second_low_onset:second_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (!is.na(first_low_onset)) {
        DS_second_low <- sum(data_year$SOLDI[first_low_onset:second_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_second_low <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= second_low_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_low_onset)) {
      DS_second_low <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY >= second_low_onset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_peak_val) && second_peak_val > low_onset_threshold) {
      DS_second_low <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_second_low <- NA
    
    # Second moderate drought severity determination
    if (exists("second_mod_offset") && !is.na(second_mod_offset)){
      if (!is.na(second_mod_onset)) {
        DS_second_mod <- sum(data_year$SOLDI[second_mod_onset:second_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (!is.na(first_mod_onset)) {
        DS_second_mod <- sum(data_year$SOLDI[first_mod_onset:second_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_second_mod <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= second_mod_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_mod_onset)) {
      DS_second_mod <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY >= second_mod_onset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_peak_val) && second_peak_val > mod_onset_threshold) {
      DS_second_mod <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_second_mod <- NA
    
    # Second extreme drought severity determination
    if (exists("second_extreme_offset") && !is.na(second_extreme_offset)){
      if (!is.na(second_extreme_onset)) {
        DS_second_extreme <- sum(data_year$SOLDI[second_extreme_onset:second_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else if (!is.na(first_extreme_onset)) {
        DS_second_extreme <- sum(data_year$SOLDI[first_extreme_onset:second_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
      } else DS_second_extreme <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY <= second_extreme_offset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_extreme_onset)) {
      DS_second_extreme <- sum(data_year$SOLDI[data_year$year == year_to_check & data_year$DOY >= second_extreme_onset]) * (data_year$DOY[2] - data_year$DOY[1])
    } else if (!is.na(second_peak_val) && second_peak_val > high_onset_threshold) {
      DS_second_extreme <- sum(data_year$SOLDI[data_year$year == year_to_check]) * (data_year$DOY[2] - data_year$DOY[1])
    } else DS_second_extreme <- NA
    
    # Reestablishing the true DOY if the drought started the previous year
    if (exists("previous_peak_doy") && !is.na(previous_peak_doy)) {
      first_peak_doy <- first_peak_doy - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_peak_doy <- second_peak_doy - length(data_year$DOY[data_year$year == year_to_check - 1])
      minDOY <- minDOY - length(data_year$DOY[data_year$year == year_to_check - 1])
      min_after_second <- min_after_second - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_low_onset <- first_low_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_low_offset <- first_low_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_mod_onset <- first_mod_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_mod_offset <- first_mod_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_extreme_onset <- first_extreme_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      first_extreme_offset <- first_extreme_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      min_B4_first <- min_B4_first - length(data_year$DOY[data_year$year == year_to_check - 1 ])
      
      second_low_onset <- second_low_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_low_offset <- second_low_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_mod_onset <- second_mod_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_mod_offset <- second_mod_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_extreme_onset <- second_extreme_onset - length(data_year$DOY[data_year$year == year_to_check - 1])
      second_extreme_offset <- second_extreme_offset - length(data_year$DOY[data_year$year == year_to_check - 1])
      
    }
    
    # Inserting the droughts facets in the final "DF" dataframe
    DF$First.Peak.Value[DF$Year== year_to_check] <- first_peak_val
    DF$First.Peak.DOY[DF$Year== year_to_check] <- first_peak_doy
    DF$First.Low.D.Onset[DF$Year== year_to_check] <- first_low_onset
    DF$First.Mod.D.Onset[DF$Year== year_to_check] <- first_mod_onset
    DF$First.Extreme.D.Onset[DF$Year== year_to_check] <- first_extreme_onset
    DF$First.Low.D.Offset[DF$Year== year_to_check] <- first_low_offset
    DF$First.Mod.D.Offset[DF$Year== year_to_check] <- first_mod_offset
    DF$First.Extreme.D.Offset[DF$Year== year_to_check] <- first_extreme_offset
    DF$First.Low.D.S[DF$Year== year_to_check] <- DS_first_low
    DF$First.Mod.D.S[DF$Year== year_to_check] <- DS_first_mod
    DF$First.Extreme.D.S[DF$Year== year_to_check] <- DS_first_extreme
    DF$Second.Low.D.S[DF$Year== year_to_check] <- DS_second_low
    DF$Second.Mod.D.S[DF$Year== year_to_check] <- DS_second_mod
    DF$Second.Extreme.D.S[DF$Year== year_to_check] <- DS_second_extreme
    DF$First.Low.D.Duration[DF$Year== year_to_check] <- first_low_duration
    DF$First.Mod.D.Duration[DF$Year== year_to_check] <- first_mod_duration
    DF$First.Extreme.D.Duration[DF$Year== year_to_check] <- first_extreme_duration
    DF$Second.Low.D.Duration[DF$Year== year_to_check] <- second_low_duration
    DF$Second.Mod.D.Duration[DF$Year== year_to_check] <- second_mod_duration
    DF$Second.Extreme.D.Duration[DF$Year== year_to_check] <- second_extreme_duration
    DF$Second.Drying.Rate[DF$Year== year_to_check] <- drying_rate
    DF$Second.Wetting.Rate[DF$Year== year_to_check] <- wetting_rate
    DF$Second.Intercept.drying[DF$Year== year_to_check] <- intercept_drying_curve
    DF$Second.Intercept.wetting[DF$Year== year_to_check] <- intercept_wetting_rate
    DF$Second.Peak.Value[DF$Year== year_to_check] <- second_peak_val
    DF$Second.Peak.DOY[DF$Year== year_to_check] <- second_peak_doy
    DF$Second.Low.D.Onset[DF$Year== year_to_check] <- second_low_onset
    DF$Second.Mod.D.Onset[DF$Year== year_to_check] <- second_mod_onset
    DF$Second.Extreme.D.Onset[DF$Year== year_to_check] <- second_extreme_onset
    DF$Second.Low.D.Offset[DF$Year== year_to_check] <- second_low_offset
    DF$Second.Mod.D.Offset[DF$Year== year_to_check] <- second_mod_offset
    DF$Second.Extreme.D.Offset[DF$Year== year_to_check] <- second_extreme_offset
    DF$Unimodal[DF$Year== year_to_check] <- unimodal
    
    # When first drought's offsets DOY and second drought's onsets DOY overlap, their values are recalculated
    if (!is.na(DF$First.Low.D.Offset[DF$Year== year_to_check ]) && !is.na(DF$Second.Low.D.Onset[DF$Year== year_to_check]) && DF$First.Low.D.Offset[DF$Year== year_to_check ] > DF$Second.Low.D.Onset[DF$Year== year_to_check]) {
      diff <- round(((abs(DF$First.Low.D.Offset[DF$Year== year_to_check] - DF$Second.Low.D.Onset[DF$Year== year_to_check]))/2) + 1)
      DF$First.Low.D.Offset[DF$Year== year_to_check] <- DF$First.Low.D.Offset[DF$Year== year_to_check] - diff
      DF$Second.Low.D.Onset[DF$Year== year_to_check] <- DF$Second.Low.D.Onset[DF$Year== year_to_check] + diff
    }
    if (!is.na(DF$First.Mod.D.Offset[DF$Year== year_to_check ]) && !is.na(DF$Second.Mod.D.Onset[DF$Year== year_to_check]) && DF$First.Mod.D.Offset[DF$Year== year_to_check ] > DF$Second.Mod.D.Onset[DF$Year== year_to_check]) {
      diff <- round(((abs(DF$First.Mod.D.Offset[DF$Year== year_to_check] - DF$Second.Mod.D.Onset[DF$Year== year_to_check]))/2) + 1)
      DF$First.Mod.D.Offset[DF$Year== year_to_check] <- DF$First.Mod.D.Offset[DF$Year== year_to_check] - diff
      DF$Second.Mod.D.Onset[DF$Year== year_to_check] <- DF$Second.Mod.D.Onset[DF$Year== year_to_check] + diff
    }    
    if (!is.na(DF$First.Extreme.D.Offset[DF$Year== year_to_check ]) && !is.na(DF$Second.Extreme.D.Onset[DF$Year== year_to_check]) && DF$First.Extreme.D.Offset[DF$Year== year_to_check ] > DF$Second.Extreme.D.Onset[DF$Year== year_to_check]) {
      diff <- round(((abs(DF$First.Extreme.D.Offset[DF$Year== year_to_check] - DF$Second.Extreme.D.Onset[DF$Year== year_to_check]))/2) + 1)
      DF$First.Extreme.D.Offset[DF$Year== year_to_check] <- DF$First.Extreme.D.Offset[DF$Year== year_to_check] - diff
      DF$Second.Extreme.D.Onset[DF$Year== year_to_check] <- DF$Second.Extreme.D.Onset[DF$Year== year_to_check] + diff
    }
    
    # Saving the last minimum value and DOY of the year necessary for next year computation
    if (!is.na(min_after_second)) {
      min_prev <- min_after_second - length(data_year$SOLDI[data_year$year == year_to_check])
      min_prev_val <- min_after_second_val
      
    } else {
      min_prev <- minDOY - length(data_year$SOLDI[data_year$year == year_to_check])
      min_prev_val <- minSOLDI 
      }
    
    # Filling the 'NA' in the final 'DF' dataframe (for onsets, offsets, and durations) in cases of multiyear droughts using values from next and previous years
    if (year_to_check == max(seq_years)) {
      # First drought onset phase
      types <- c("Low", "Mod", "Extreme")
      thresholds <- c(low_onset_threshold, mod_onset_threshold, high_onset_threshold)
      
      for (i in seq_along(types)) {
        type <- types[i]
        threshold <- thresholds[i]
        
        first_col <- paste0("First.", type, ".D.Onset")
        second_col <- paste0("Second.", type, ".D.Onset")
        
        for (n in 2:nrow(DF)) {
          if (is.na(DF[[first_col]][n]) && !is.na(DF$First.Peak.Value[n]) && DF$First.Peak.Value[n] > threshold) {
            idx <- n - 1
            
            while (idx > 0) {
              if (!is.na(DF[[second_col]][idx])) {
                DF[[first_col]][n] <- DF[[second_col]][idx] - 365
                break
              } else if (!is.na(DF[[first_col]][idx])) {
                DF[[first_col]][n] <- DF[[first_col]][idx] - 365
                break
              }
              idx <- idx - 1
            }
          }
        }
      }
      
      # Second drought onset phase      
      for (i in seq_along(types)) {
        type <- types[i]
        threshold <- thresholds[i]
        
        first_onset <- paste0("First.", type, ".D.Onset")
        second_onset <- paste0("Second.", type, ".D.Onset")
        second_peak <- "Second.Peak.Value"
        
        for (n in 1:nrow(DF)) {
          if (is.na(DF[[second_onset]][n]) && !is.na(DF[[second_peak]][n]) && DF[[second_peak]][n] > threshold) {
            
            if (!is.na(DF[[first_onset]][n])) {
              DF[[second_onset]][n] <- DF[[first_onset]][n]
              
            } else {
              idx <- n - 1
              
              while (idx > 0) {
                if (!is.na(DF[[second_onset]][idx])) {
                  DF[[second_onset]][n] <- DF[[second_onset]][idx] - 365
                  break
                } else if (!is.na(DF[[first_onset]][idx])) {
                  DF[[second_onset]][n] <- DF[[first_onset]][idx] - 365
                  break
                }
                idx <- idx - 1
              }
            }
          }
        }
      }
      
      # First drought offset phase
      thresholds <- c(low_offset_threshold, mod_offset_threshold, high_offset_threshold)
      
      for (i in seq_along(types)) {
        type <- types[i]
        threshold <- thresholds[i]
        
        first_offset <- paste0("First.", type, ".D.Offset")
        second_offset <- paste0("Second.", type, ".D.Offset")
        peak_SOLDI <- "First.Peak.Value"

        for (n in 1:(nrow(DF))) {
          if (is.na(DF[[first_offset]][n]) && !is.na(DF[[peak_SOLDI]][n]) && DF[[peak_SOLDI]][n] > threshold) { # !!!! ici il faut aussi verifier que le minimum est en dssous de threshold sinon rerecherche de offset 
            
            if (!is.na(DF[[second_offset]][n])) {
              DF[[first_offset]][n] <- DF[[second_offset]][n]
              
            } else {
              idx <- n + 1
              
              found <- FALSE
              while (idx <= nrow(DF)) {
                if (!is.na(DF[[first_offset]][idx])) {
                  DF[[first_offset]][n] <- DF[[first_offset]][idx] + 365 * (idx - n)
                  found <- TRUE
                  break
                } else if (!is.na(DF[[second_offset]][idx])) {
                  DF[[first_offset]][n] <- DF[[second_offset]][idx] + 365 * (idx - n )
                  found <- TRUE
                  break
                }
                
                idx <- idx + 1
              }
              
              if (!found) {
                DF[[first_offset]][n] <- 365 * (idx- n)
              }
              
            }
          }
        }
      }

      # Second drought offset phase
      for (i in seq_along(types)) {
        type <- types[i]
        threshold <- thresholds[i]
        
        first_offset <- paste0("First.", type, ".D.Offset")
        second_offset <- paste0("Second.", type, ".D.Offset")
        peak_SOLDI <- "Second.Peak.Value"
        
        for (n in 1:(nrow(DF) - 1)) {
          if (is.na(DF[[second_offset]][n]) && !is.na(DF[[peak_SOLDI]][n]) && DF[[peak_SOLDI]][n] > threshold && DF$Unimodal[n] == FALSE) {
            idx <- n + 1
            
            while (idx <= nrow(DF)) {
              if (!is.na(DF[[first_offset]][idx])) {
                DF[[second_offset]][n] <- DF[[first_offset]][idx] + 365 * (idx - n)
                break
              } else if (!is.na(DF[[second_offset]][idx])) {
                DF[[second_offset]][n] <- DF[[second_offset]][idx] + 365 * (idx - n)
                break
              } else DF[[second_offset]][n] <- 365
              idx <- idx + 1
            }
          }
        }
      }

      # New drought duration calculations
      for (n in 1:nrow(DF)) {
        if (is.na(DF$First.Low.D.Duration[n]) && !is.na(DF$First.Peak.Value[n]) && DF$First.Peak.Value[n] > low_onset_threshold) {
          DF$First.Low.D.Duration[n] <- DF$First.Low.D.Offset[n] - DF$First.Low.D.Onset[n] + 1
        }
        if (is.na(DF$First.Mod.D.Duration[n]) && !is.na(DF$First.Peak.Value[n]) && DF$First.Peak.Value[n] > mod_onset_threshold) {
          DF$First.Mod.D.Duration[n] <-  DF$First.Mod.D.Offset[n] - DF$First.Mod.D.Onset[n] + 1
        }
        if (is.na(DF$First.Extreme.D.Duration[n]) && !is.na(DF$First.Peak.Value[n]) && DF$First.Peak.Value[n] > high_onset_threshold) {
          DF$First.Extreme.D.Duration[n] <-  DF$First.Extreme.D.Offset[n] - DF$First.Extreme.D.Onset[n] + 1
        }
        
        if (is.na(DF$Second.Low.D.Duration[n]) && !is.na(DF$Second.Peak.Value[n]) && DF$Second.Peak.Value[n] > low_onset_threshold) {
          DF$Second.Low.D.Duration[n] <-  DF$Second.Low.D.Offset[n] - DF$Second.Low.D.Onset[n] + 1
        }
        if (is.na(DF$Second.Mod.D.Duration[n]) && !is.na(DF$Second.Peak.Value[n]) && DF$Second.Peak.Value[n] > mod_onset_threshold) {
          DF$Second.Mod.D.Duration[n] <-  DF$Second.Mod.D.Offset[n] - DF$Second.Mod.D.Onset[n] + 1
        }
        if (is.na(DF$Second.Extreme.D.Duration[n]) && !is.na(DF$Second.Peak.Value[n]) && DF$Second.Peak.Value[n] > high_onset_threshold) {
          DF$Second.Extreme.D.Duration[n] <-  DF$Second.Extreme.D.Offset[n] - DF$Second.Extreme.D.Onset[n] + 1
        }
        
      # Monsoon duration calculations
        if (!is.na(DF$First.Low.D.Offset[n])){
          if (!is.na(DF$Second.Low.D.Onset[n]) && DF$Second.Low.D.Onset[n] > DF$First.Low.D.Offset[n]) {
            DF$Wet.Monsoon.Duration[n] <- DF$Second.Low.D.Onset[n] - DF$First.Low.D.Offset[n] +1
          } else if (!is.na(DF$Second.Peak.DOY[n]) && DF$Second.Peak.DOY[n] > DF$First.Low.D.Offset[n]) {
            DF$Wet.Monsoon.Duration[n] <- DF$Second.Peak.DOY[n] - DF$First.Low.D.Offset[n]+1
          } else if (!is.na(DF$Second.Peak.DOY[n])){ 
            DF$Wet.Monsoon.Duration[n] <- 0
            } 
        } else if (!is.na(DF$Second.Low.D.Onset[n]) && DF$Second.Low.D.Onset[n] > DF$First.Peak.DOY[n]) {
          DF$Wet.Monsoon.Duration[n] <- DF$Second.Low.D.Onset[n] - DF$First.Peak.DOY[n]+1
        } else if (!is.na(DF$Second.Extreme.D.Onset[n])){
          DF$Wet.Monsoon.Duration[n] <- 0 
        }
        
        if (!is.na(DF$First.Mod.D.Offset[n])){
          if (!is.na(DF$Second.Mod.D.Onset[n]) && DF$Second.Mod.D.Onset[n] > DF$First.Mod.D.Offset[n]) {
            DF$Mod.Monsoon.Duration[n] <- DF$Second.Mod.D.Onset[n] - DF$First.Mod.D.Offset[n] +1
          } else if (!is.na(DF$Second.Peak.DOY[n]) && DF$Second.Peak.DOY[n] > DF$First.Mod.D.Offset[n]) {
            DF$Mod.Monsoon.Duration[n] <- DF$Second.Peak.DOY[n] - DF$First.Mod.D.Offset[n]+1
          } else if (!is.na(DF$Second.Peak.DOY[n])){ 
            DF$Mod.Monsoon.Duration[n] <- 0
          } 
        } else if (!is.na(DF$Second.Mod.D.Onset[n]) && DF$Second.Mod.D.Onset[n] > DF$First.Peak.DOY[n]) {
          DF$Mod.Monsoon.Duration[n] <- DF$Second.Mod.D.Onset[n] - DF$First.Peak.DOY[n]+1
        } else if (!is.na(DF$Second.Extreme.D.Onset[n])){
          DF$Mod.Monsoon.Duration[n] <- 0 
        }
        
        if (!is.na(DF$First.Extreme.D.Offset[n])){
          if (!is.na(DF$Second.Extreme.D.Onset[n]) && DF$Second.Extreme.D.Onset[n] > DF$First.Extreme.D.Offset[n]) {
            DF$Dry.Monsoon.Duration[n] <- DF$Second.Extreme.D.Onset[n] - DF$First.Extreme.D.Offset[n] +1
          } else if (!is.na(DF$Second.Peak.DOY[n]) && DF$Second.Peak.DOY[n] > DF$First.Extreme.D.Offset[n]) {
            DF$Dry.Monsoon.Duration[n] <- DF$Second.Peak.DOY[n] - DF$First.Extreme.D.Offset[n]+1
          } else if (!is.na(DF$Second.Peak.DOY[n])){ 
            DF$Dry.Monsoon.Duration[n] <- 0
          } 
        } else if (!is.na(DF$Second.Extreme.D.Onset[n]) && DF$Second.Extreme.D.Onset[n] > DF$First.Peak.DOY[n]) {
          DF$Dry.Monsoon.Duration[n] <- DF$Second.Extreme.D.Onset[n] - DF$First.Peak.DOY[n]+1
        } else if (!is.na(DF$Second.Extreme.D.Onset[n])){
          DF$Dry.Monsoon.Duration[n] <- 0 
          }
      
      }
    }
  }
  
  # Export the final dataframe of the pixel
  res <- tryCatch({
    write.csv(DF, file = paste0(output_dir,"DF_DRY_NAM_", coord_x_y, ".csv"), row.names = FALSE)
    TRUE
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    FALSE
  })
  return(res)
}



#### 3.a. Execute the functions with 'lapply' ####

# lapply(files_to_process, function(f) {
#   process_file(f, date, seq_years, DF, low_onset_threshold, low_offset_threshold,
#                mod_onset_threshold,mod_offset_threshold,
#                high_onset_threshold, high_offset_threshold, output_dir)
# })

#### 3.b. Execute the functions with 'parallel' ####

nb_cores <- detectCores() - 10
cl <- makeCluster(nb_cores)

# Export function and variables
clusterExport(cl, varlist = c("process_file","peakfunction", "date", "seq_years", "DF", "low_onset_threshold", "low_offset_threshold",
                              "mod_onset_threshold","mod_offset_threshold",
                              "high_onset_threshold", "high_offset_threshold", "output_dir"))

# Load necessary packages
clusterEvalQ(cl, {
  library(terra)
  library(lubridate)
  library(dplyr)
  library(tidyr)
  library('phenofit')
})

# Run the parallel process
results_parallel <- parLapply(cl, files_to_process, function(f) {
  process_file(f, date, seq_years, DF, low_onset_threshold, low_offset_threshold,
               mod_onset_threshold,mod_offset_threshold,
               high_onset_threshold, high_offset_threshold, output_dir)
})


stopCluster(cl)





















