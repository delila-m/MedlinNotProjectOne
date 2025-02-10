library(raster)
library(sf)
library(terra)
library(tigris)
library(dplyr)
library(tidyverse)
library(tidyterra)
library(ncdf4)
library(randomForest)
library(caret)

# read in the data, do some cleaning.
# First the usdm data for coconino county

drought.az.cocoinino <- read.csv("data-raw/USDM-04005.csv")


# now the pdsi for the whole US
# pdsi data
pdsi <- rast("data-raw/agg_met_pdsi_1979_CurrentYear_CONUS.nc")

# This function crops the pdsi file to the county level
# Takes the state and county desired as strings, and a spatraster object containing pdsi data
## note this raster file has to be loaded in using the rast function from the terra package
crop.county.pdsi <- function(state, county, spat.us.pdsi){ # county and state need to be saved as strings
  # get the specific county
  sf.county <- counties(state = state,
                        class = "sf")
  # first I'll start with just selecting  county
  county <- sf.county[sf.county$NAME == county, ]

  # use extract function, cells and xy argument return the cell number and xy coordinates, respectively.
  # added to help with mapping visualizatona later on
  pdsi.extracted <- terra::extract(spat.us.pdsi, county,
                                   cells = TRUE, xy = TRUE) # creates a data frame

  return(pdsi.extracted)
}

# This function normalizes the land area percentages for the US drought index
# note this should be done before calculating the weighted averages
normalize.usdm <- function(data){
  # D4 stays the same
  # D3 gets D4 subtracted from it
  data$normal.D3 <- data$D3 - data$D4
  # D2 has D3 subtracted from it
  data$normal.D2 <- data$D2 - data$D3
  # and so on
  data$normal.D1 <- data$D1 - data$D2
  data$normal.D0 <- data$D0 - data$D1
  # save values from temporary columns to the original to keep names consistent
  data$D3 <- data$normal.D3
  data$D2 <- data$normal.D2
  data$D1 <- data$normal.D1
  data$D0 <- data$normal.D0

  # return the data frame minus the temporary columns
  data <- subset(data,
                 select = -c(normal.D3, normal.D2, normal.D1, normal.D0))
  # should we test leaving the ugly drought percentages in?
  return(data)
}

# This function calculates and creates a new column for the weighted average of the US drought index
# the average is weighted by the land area percent for each index
usdm.weighted.average <- function(data) {
  # Define the weights
  weights <- 0:4

  # Calculate the weighted average across columns D0 to D4
  data$USDM_Avg <- mapply(function(D0, D1, D2, D3, D4) {
    values <- c(D0, D1, D2, D3, D4)
    weighted_sum <- sum(values * weights)
    weight_sum <- sum(values)

    if (weight_sum == 0) {
      return(NA)  # Avoid division by zero error
    } else {
      return(weighted_sum / weight_sum)
    }
  }, data$D0, data$D1, data$D2, data$D3, data$D4)

  return(data)  # Return the modified data frame
}

# This function interpolates values for PDSI to find the midpoint of the 5 day cycle which the data is collected on.
# the midpoint date is the same for the interpolate.usdm function
interpolate.pdsi <- function(pdsi.data){

  # assign grid cell IDs
  pdsi.data$ID <- seq_len(nrow(pdsi.data))

  # select only the layers which have pdsi measurements
  goodData <- select(pdsi.data, ID, cell, x, y,
                     starts_with("daily_mean_palmer_drought_severity_index_day="))

  # get rid of of the long name
  allNames <- names(goodData)
  newNames <- str_remove(allNames,
                         "daily_mean_palmer_drought_severity_index_day=")
  names(goodData) <- newNames

  # pivot so we can have all of the dates in one column
  goodLong <- pivot_longer(goodData, cols = -c(ID, x, y, cell),
                           names_to = "Day", values_to = "PDSI") %>%
    # deal with the date
    mutate(
      day = as.integer(Day),
      Date = as.Date(day, origin = "1900-01-01"), # starting date
      midPointDay = day - 2.5) # subtract 2.5 to find the midpoint of the 5 day cycle

  return(goodLong)
}

# This function interpolates values for the us drought monitor
# the interpolated dates are the same for the PDSI data
interpolate.usdm <- function(drought.data){
  # deal with date
  baselineDate <- as.Date("1900-01-01") # still didn't find the start date from NOAA..

  # now we need to convert the drought monitor into days since 1900, and calculate midPointDay
  drought.data$ValidMid <- as.Date(drought.data$ValidEnd) - 3.5 # subtract 3.5 to find the midpoint of the 7 day cycle

  # actually apply the dates to our drought data
  drought.data$midPointDay <- as.numeric(drought.data$ValidMid - baselineDate)

  return(drought.data)
}

# This function joins the interpolated USDM and PDSI data sets, and puts the data in a format ready for modeling
join.pdsi.usdm <- function(pdsi.data, usdm.data){
  # we simply have to join by date
  drought.pdsi <- left_join(pdsi.data, usdm.data)

  # leaving out any na values
  drought.pdsi <- na.omit(drought.pdsi)

  # select only the columns we want
  reshaped_data <- drought.pdsi %>%
    select(Date, x, y, cell, USDM_Avg, PDSI)

  return(reshaped_data)
}

# This function combines all of the above cleaning for a streamlined approach
# the arguments needed are the state and county names as strings,
# the pdsi dataset from drought.gov, loaded in as a spatraster object,
# and the usdm data from drought.gov for the specific county
clean.county.data <- function(state, county, spat.us.pdsi, cropped.usdm.data){
  # crop the pdsi data
  pdsi <- crop.county.pdsi(state, county, spat.us.pdsi)

  # clean the usdm data
  normal.drought <- normalize.usdm(cropped.usdm.data)
  weighted.drought <- usdm.weighted.average(normal.drought)
  interpolated.drought <- interpolate.usdm(weighted.drought)

  # clean the pdsi data
  interpolated.pdsi <- interpolate.pdsi(pdsi)

  # join the datasets and reshaping
  full.data <- join.pdsi.usdm(interpolated.pdsi, interpolated.drought)

  # return the final cleaned dataset
  return(full.data)
}


# Use all of the above functions to clean the data
Coconino_USDM_PDSI <- clean.county.data("Arizona", "Coconino", pdsi, drought.az.cocoinino)
