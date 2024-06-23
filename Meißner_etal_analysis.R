#####################    WildLive! Jaguar Proj.  #####################
##################### Mapping & Prelim. Analysis #####################

## Install packages first! 
# remotes::install_github("carlopacioni/camtrapRdeluxe")

library(camtrapR)   # Camera traps
library(secr)       # SECR Analysis
library(tidyverse)  # Graphics and data manipulation
library(activity)   # Activity curves 
library(data.table) # Data manipulation
library(lubridate)  # Dates 
library(ggmap)      # Base maps
library(ggspatial)  # Mapping
library(ggsn)       # Mapping
library(scatterpie) # Mapping
library(geosphere)  # Mapping
library(scales)     # Mapping
library(gridExtra)  # Mapping 
library(GISTools)   # Mapping
library(sp)         # Mapping
library(sf)         # Mapping
library(ggrepel)    # Labels
library(mapview)    # Interactive display
library(pryr)       # Computations
library(rstudioapi) # API
library(viridis)    # Aesthetics
library(ggpubr)     # Plot arrangements 
library(remotes)    # misc
library(gdata)      # misc
library(miscTools)  # misc
library(extrafont)  # fonts
library(cowplot)    # plot grid
library(Hmisc)      # misc
library(rstatix)    # hypothesis tests
library(datarium)   # hypothesis tests 

setwd("") # set WD 
list.files() # check data availability
options(scipen = 999) # disable sci. notation

################# (1) Data cleaning & preparation #################

## ----- 1.1. Read & clean record data -----

## Load Metadata & set formats 
dat.all <- read.csv("jaguar_database.csv", header = TRUE, sep = ",", na.strings = c("", "NA")) # import data

## Clean data
dat.clean <- dat.all[, c(6,15,12,13,24,25,26,27,23,22)]   # relevant columns 
colnames(dat.clean) <- c("ID", "DateTimeOriginal", "Age", "Sex", "Station", "Camera", "lat", "long","Grid", "Setup")

setDT(dat.clean) # check and remove missing data 
missing  <- na.omit(dat.clean, cols = c("DateTimeOriginal", "Station"), invert = TRUE)
nrow(subset(missing, is.na(missing$Station)==TRUE)) # 1 unidentified station
nrow(subset(missing, is.na(missing$DateTimeOriginal)==TRUE)) # 127 captures lack time stamp

## Remove values with missing data in relevant columns 
dat.clean <- na.omit(dat.clean, cols = c("DateTimeOriginal", "Station"), invert = FALSE)
nrow(dat.all)-nrow(dat.clean) #127 rows lack DateTimeOriginal 

## Make sure R understands DateTimeOriginal is a date format
DateTimeOriginal <- as.POSIXct(dat.clean$DateTimeOriginal,
                               format = "%Y:%m:%d %H:%M:%S")

## Remove all potential white spaces in data 
dat.clean <- as.data.frame(apply(dat.clean,2,function(x)gsub('\\s+', '',x)))
dat.clean$DateTimeOriginal <- DateTimeOriginal # re-add DateTimeOriginal

# We exclude all observations from 2020 onwards due to technical problems. 
dat.clean <- subset(dat.clean, DateTimeOriginal <= "2020-01-01 00:00:00")

## ----- 1.2. Check coordinates & area size -----

## Build a spatial object to make a preliminary map 
dat.clean$lat  <- as.numeric(dat.clean$lat) # make sure formats are correct
dat.clean$long <- as.numeric(dat.clean$long) # make sure formats are correct
coords <- st_as_sf(dat.clean[ ,c("Station", "lat", "long")], # build spatial object
                   coords = c("long", "lat"), crs = 4326)    # using WGS84 (4326)
mapview(coords, legend = TRUE) # examine points on map 


############ (2) Activity curves & analysis ############ 

# For the activity curves we use only captures of adult individuals 
# which are temporally independent by a delta time of 5 hours (300 min). 

require(camtrapRdeluxe, quietly = TRUE, warn.conflicts = FALSE) # get delta time 
actRTI <- assessTemporalIndependence(intable = dat.clean,      
                                                    deltaTimeComparedTo = "lastIndependentRecord", 
                                                    columnOfInterest    = "ID",
                                                    stationCol          = "Station", 
                                                    camerasIndependent  = FALSE, 
                                                    minDeltaTime        = 60)


## Next we convert the time stamps in radian time (on the range [0, 2*pi])
act_time <- gettime(actRTI$DateTimeOriginal, "%Y-%m-%d %H:%M:%S")

## ----- 2.1. All relevant individuals -----

# We can fit an activity model and get meaningful confidence intervals 
# by bootstrapping. This can either be done by sampling the data or the 
# fitted distribution. The best selection often depends on the available sample size. 
# If the sample size is sufficient, we can sample from the data.

jag_act <- fitact(act_time, reps = 100, sample = "data")

## Examine activity model 
plot(jag_act, yunit="frequency", data = "both",
     tline=list(col="red", lwd = 2), 
     cline=list(col="red"))

## For aesthetic purposes... plot the curve in ggplot2. 
activity <- data.frame(jag_act@pdf)

## Compute y values to "frequency"
activity$se <- NULL # drop standard error 
activity$y <- activity$y*length(jag_act@data)*pi/12     # frequency
activity$lcl <- activity$lcl*length(jag_act@data)*pi/12 # lower confidence interval
activity$ucl <- activity$ucl*length(jag_act@data)*pi/12 # upper confidence interval
activity$x <- activity$x*12/pi                          # time 
colnames(activity) <- c("time", "frequency", "lcl","ucl")

hist <- data.frame(jag_act@data*12/pi)
colnames(hist) <- "time"
hist$time <- as.numeric(hist$time)

act_plot <- ggplot()+
  # Plot fit model 
  geom_line(data = activity, 
            aes(x=time, y=frequency, 
                color = "Activity", linetype = "Activity"), size = .3) +
  # Upper confidence interval 
  geom_line(data = activity, 
            aes(x=time, y=ucl, color = "CI", 
                linetype = "CI"), size = .3) +
  # Lower confidence interval
  geom_line(data = activity, 
            aes(x=time, y=lcl, color = "CI", 
                linetype = "CI"), size = .3) +
  # Raw data histogram
  geom_histogram(data = hist, aes(x=time, fill = "Raw data"),
                 color = "transparent", alpha = .15, binwidth = 1) +
  # Set color & fill aesthetics 
  scale_fill_manual(name = "", values = c("Raw data" = "grey20")) +
  scale_color_manual(name = "", values = c("Activity" = "black",
                                           "CI" = "grey50")) +
  scale_linetype_manual(name = "", values = c("Activity" = "solid",
                                              "CI" = "dotted")) +
  # X-axis breaks 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4))+
  # Theme aesthetics 
  theme_bw(base_line_size = .2)+
  theme(axis.title.x=element_text(family="Arial", size = 6),
        axis.title.y=element_text(family="Arial", size = 6), 
        axis.text.x=element_text(family="Arial", size = 6),
        axis.text.y=element_text(family="Arial", size = 6),
        # Legend aesthetics 
        legend.text=element_text(family="Arial", size = 6),
        legend.key.width = unit(.25, 'cm'),
        legend.key.height = unit(.05, 'cm'),
        legend.background = element_rect(fill="white", size=0.15, linetype="solid", colour ="black"),
        legend.margin = margin(t=.175, r=.125, b=0, l=.125, unit = "cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.position = c(.7325,.915),
        legend.spacing.x = unit(0.125, "lines"))+
  # Labels 
  ylab("Frequency") +
  xlab("Time of day")

# Save & export activity plot 
ggsave(filename = "fig_3.png",
       plot = act_plot,
       bg = "white",
       width = 80,
       height = 50,
       units = "mm",
       dpi = 1257)

## ------ 2.2 Activity curves for individuals separated ------

# We also check individuals separately if they have sufficient captures

table(actRTI$ID) # only 1,2,4 and 5 make sense.

jag_1 <- subset(actRTI, ID == "F-01") # get individual data
jag_2 <- subset(actRTI, ID == "M-01") # get individual data
jag_4 <- subset(actRTI, ID == "F-03") # get individual data
jag_5 <- subset(actRTI, ID == "M-02") # get individual data

jag_1 <- gettime(jag_1$DateTimeOriginal, "%Y-%m-%d %H:%M:%S") # get radian time
jag_2 <- gettime(jag_2$DateTimeOriginal, "%Y-%m-%d %H:%M:%S") # get radian time
jag_4 <- gettime(jag_4$DateTimeOriginal, "%Y-%m-%d %H:%M:%S") # get radian time
jag_5 <- gettime(jag_5$DateTimeOriginal, "%Y-%m-%d %H:%M:%S") # get radian time

jag_1 <- fitact(jag_1, reps = 100, sample = "model") # bootstrap by fitted values
jag_2 <- fitact(jag_2, reps = 100, sample = "model") # bootstrap by fitted values
jag_4 <- fitact(jag_4, reps = 100, sample = "model") # bootstrap by fitted values
jag_5 <- fitact(jag_5, reps = 100, sample = "model") # bootstrap by fitted values

par(mfrow = c(2, 2)) # examine curves 
plot(jag_1, yunit="frequency", data = "both", tline=list(col="red", lwd = 2), cline=list(col="red"), main = "F-01")
plot(jag_2, yunit="frequency", data = "both", tline=list(col="red", lwd = 2), cline=list(col="red"), main = "M-01")
plot(jag_4, yunit="frequency", data = "both", tline=list(col="red", lwd = 2), cline=list(col="red"), main = "F-03")
plot(jag_5, yunit="frequency", data = "both", tline=list(col="red", lwd = 2), cline=list(col="red"), main = "M-02")

# Clean Global Environment 
rm(jag_1, jag_2, jag_4, jag_5, activity, actRTI, coords, missing, DateTimeOriginal, hist, jag_act, act_time)


############ (3) Spatially explicit capture recapture analysis ############ 

## ------ 3.1 Capture data ------

# For the SECR analysis we use only captures of (sub)adults which 
# are temporally independent by at least 60 minutes 

recordTableIndividual <- assessTemporalIndependence(intable = subset(dat.clean,
                                                                     ID != "J-01-1" &  
                                                                     ID != "J-03-1" &
                                                                     ID != "J-03-3"),      
                                                    deltaTimeComparedTo = "lastIndependentRecord", 
                                                    columnOfInterest    = "ID",
                                                    stationCol          = "Station", 
                                                    camerasIndependent  = FALSE, 
                                                    minDeltaTime        = 60)
detach("package:camtrapRdeluxe", unload = TRUE)
recordTableIndividual$Station <- str_trim(recordTableIndividual$Station, side = "both")
recordTableIndividual$Species = "Jaguar"


## ------ 3.2 Camera trap data ------

# We need to clean the camtrap data and make sure that the relevant variables
# strictly match in the capture and camtrap data files.

# Import camtrap data 
camtraps <- read.csv("camtraps.csv", header = TRUE, na.strings = "", sep = ",") 

# Remove all potential white spaces in data 
camtraps <- as.data.frame(apply(camtraps, 2, function(x)gsub('\\s+', '',x)))
camtraps$Station <- str_trim(camtraps$Station, side ="both")

## Check if Station, Camera and Session in 'camtraps' match 'recordTableIndividual'
match_cam <- camtraps %>% dplyr::select(Station) %>% arrange() %>% distinct() 
match_rec <- recordTableIndividual %>% dplyr::select(Station) %>% arrange() %>% distinct()

setDT(match_cam)  # set data table 
setDT(match_rec)  # set data table 

shared <- subset(match_cam, Station %in% match_rec$Station)     # camtraps in recordtable
missing <- subset(match_cam, !(Station %in% match_rec$Station)) # camtraps not in recordtable

# no observations from stations...
paste(missing$Station)
# ...in recordTableIndividual. 

setDT(camtraps) # remove missing values (should be none)
camtraps <- na.omit(camtraps, cols = c("Station","Setup_date"))

# ~~~~~ From here onward we keep only stations in new grid ~~~~~

## ( !! this affects the RecordTableIndividual data as well !! )

# Get only camtraps within the grid (thus stations where NewGrid == y)
grid_traps <- subset(camtraps, NewGrid == "y")

## Double check that there are ABSOLUTELY NO F***?NG white spaces in the relevant columns 
grid_traps$Station <- str_trim(grid_traps$Station, side = "both")
recordTableIndividual$Station <- str_trim(recordTableIndividual$Station, side = "both")

## Make sure only Stations present in both data frames are used !!
grid_traps <- subset(grid_traps, Station %in% recordTableIndividual$Station) # keep only Stations found in record data
recordTableIndividual <- subset(recordTableIndividual, Station %in% grid_traps$Station) # keep only Stations in grid 

grid_traps$Setup_date <- as.Date(grid_traps$Setup_date, format = "%d.%m.%y")         # set formats
grid_traps$Retrieval_date <- as.Date(grid_traps$Retrieval_date, format = "%d.%m.%y") # set formats

## Generate Camera Operation data for grid setup 
camop_problem <- camtrapR::cameraOperation(CTtable      = grid_traps,
                                           stationCol   = "Station",
                                           setupCol     = "Setup_date",
                                           cameraCol    = "Camera",
                                           retrievalCol = "Retrieval_date",
                                           writecsv     = FALSE,
                                           hasProblems  = FALSE,
                                           byCamera     = FALSE,
                                           allCamsOn    = FALSE,
                                           camerasIndependent = FALSE,
                                           dateFormat   = "%Y-%m-%d")

## Disregard multiple cameras per stations for building the capthist object
grid_traps <- distinct(grid_traps, Station, .keep_all = TRUE) 

## Set grid camtraps spatial projection
grid_traps$lat <- as.numeric(grid_traps$lat)   # set format
grid_traps$long <- as.numeric(grid_traps$long) # set format 
grid_coords <- st_as_sf(grid_traps[ ,c("Station", "lat", "long")],
                        coords = c("long", "lat"), crs = 4326) 
mapview(grid_coords) # check with map view

grid_traps            <- as.data.frame(grid_traps)            # set df
recordTableIndividual <- as.data.frame(recordTableIndividual) # set df 

## Quickly check the area size of the "new grid"" setup
areaPolygon(grid_traps[,c("lat","long")])*0.0001 # 970.4593 ha

## To build a capthist object via spatialDetectionHistory() 
## we require projected coordinates in UTM.

# First we make sure R understands that the current grid coordinates are in lat long 
grid.utm <- SpatialPoints(cbind(grid_traps$long, grid_traps$lat), 
                          proj4string = CRS("+proj=longlat"))

# Then we set the UTM Bolivia projection
grid.utm <- spTransform(grid.utm, CRS("+init=epsg:32719"))
grid.utm <- as.data.frame(grid.utm@coords) # get coordinates
colnames(grid.utm) <- c("x", "y") # set column names 

# Add UTM coordinates to grid data
grid_traps$x <- grid.utm$x
grid_traps$y <- grid.utm$y


## ---------- 3.3 Check if all variables match ----------

## First we check if the oordinates match in 
## recordTableIndividual (RC) and the camtrap table (CT)

CT <- dplyr::select(grid_traps, Station, lat, long)
RC <- dplyr::select(recordTableIndividual, Station, lat, long)

RC <- RC %>% distinct(Station, .keep_all = TRUE) %>% arrange(Station)
CT <- CT %>% distinct(Station, .keep_all = TRUE) %>% arrange(Station)

identical(RC$Station, CT$Station) # TRUE
identical(RC$long, CT$long)       # TRUE
identical(RC$lat, CT$lat)         # TRUE

# If all TRUE, they match. 

## Next, we check if the names of the stations match

CT <- grid_traps %>% dplyr::select(Station) %>% distinct(Station) %>% arrange(Station)
RC <- recordTableIndividual %>% dplyr::select(Station) %>% distinct(Station) %>% arrange(Station)
OP <- data.frame(row.names(camop_problem))
colnames(OP) <- "Station"
OP <- arrange(OP, Station)

identical(CT$Station, OP$Station) # TRUE
identical(RC$Station, OP$Station) # TRUE 

## If all TRUE, they match. 

## Lastly, we set the desired data frames and do one last check. 

grid_traps <- as.data.frame(grid_traps)    # set df
grid_traps <- arrange(grid_traps, Station) # set df

identical(rownames(camop_problem), unique(grid_traps$Station)) # TRUE

## Clean global environment
rm(grid_coords, missing, OP, RC, CT, match_cam, match_rec, shared, grid.utm)
dev.off()

## --------------- 3.4 Identifying slots for estimating jaguar density --------------------

## As we need to assume a closed population during the time to which we apply our SECR procedure, 
## we need to isolate suitable slots of ~180 days (see Harmsen et al., 2020, https://doi.org/10.1317/journal.pone.0227468)
## which feature sufficient captures of jaguars at as many stations possible. 

## Harmsen et al. argue that 180 days in comparison with traditionally used 90 days are more stable 
## and more precise in long-lived animals. (See also https://doi.org/10.1111/2041-210X.13158)

# Taking grid camera activity into account we only have 24.02.2019 - 31.12.2020 at our disposal
plot_grid( ggplot(subset(recordTableIndividual,
                         DateTimeOriginal >= "2019-01-01 00:00:00" & DateTimeOriginal <= "2019-12-31 23:59:59"))+
             geom_point(aes(x=DateTimeOriginal, y=ID))+
             theme_gray()+
             xlab("")+
             ylab(""),
           
           ggplot(subset(recordTableIndividual,
                         DateTimeOriginal >= "2019-01-01 00:00:00" & DateTimeOriginal <= "2019-12-31 23:59:59"))+
             geom_point(aes(x=DateTimeOriginal, y=Station))+
             theme_gray()+
             xlab("")+
             ylab(""),
           
           labels = c("Individuals", "Stations")
)

## ------------------- 3.5 SECR Analysis per slot  --------------------

## ~~~~~~~~~~~~~~~~~~~~~ SLOT A ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_a <- subset(recordTableIndividual, DateTimeOriginal >= "2019-01-01 00:00:01" & DateTimeOriginal <= "2019-06-29 23:59:59")

## Generate capthist object while considering sampling effort 
CT_a <- subset(grid_traps, Station %in% slot_a$Station) # get CT table for slot A
camop_a <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_a))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_a))` >= "2019-01-01 00:00:01" & `as.Date(colnames(camop_a))` <= "2019-06-29 23:59:59")
dates <- colnames(camop_a) # get desired dates 
camop_a <- as.data.frame(t(camop_a)) # prepare matrix 
camop_a$date <- rownames(camop_a) # add dates as column for sub-setting
camop_a <- subset(camop_a, date %in% as.character(slot$`as.Date(colnames(camop_a))`)) # keep only dates in slot A
camop_a[,"date"] <- NULL # drop redundant dates
camop_a <- as.matrix(t(camop_a)) # rearrange as matrix
camop_a <- camop_a[rownames(camop_a) %in% slot_a$Station, ] # subset by stations found in record table 
slot_a$Species <- "jaguar"
# Build capthist object
cap_a <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_a,
                                           CTtable = CT_a,
                                           camOp = camop_a,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-01-01",
                                           includeEffort = TRUE)

## Examine capture history and trapping effort

summary(cap_a) # examine capthist
usage(traps(cap_a)) # examine trapping effort per occasion

plot(cap_a, tracks = TRUE, main = "Slot A") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_a))
hist(m, xlab = "Movement m", main = "Slot A")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_a, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_a, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_a), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_a <- secr.fit(cap_a, buffer = 5*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_a # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_a, limits = FALSE) # plot detection probability

esa.plot(fit_a) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') # increase buffer slightly 
abline(v = 6 * initialsigma, lty = 2, col = 'green') # better

# Adjust mask buffer size 
mask <- make.mask(traps(cap_a), buffer = 5*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_a.HN <- secr.fit(cap_a, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_a.EX <- secr.fit(cap_a, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_a.HN, EX = fit_a.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 6 * initialsigma)

# While EX provides a lower AIC, the D estimates of both functions are very similar.
# HN should probably be preferred as it appears less sensitive to buffer width. 

# While we could check for different model structures, the SECR developers states, that
# "the returns from exhaustively pursuing the last silver improvement in fit are usually
# trivial". The default SECR model (i.e., D, g0 and sigma are constant) should be robust enough.

# Therefore I think using a default model structure is sufficient (as running different model  
# structures also takes a lot of processing time). 

# Also, I do not think we need to address learned responses, i.e., that the experience of 
# a capture affects the capture probability the next capture. 

# This should be the final model for slot A:

fit_a <- secr.fit(cap_a, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")


## ~~~~~~~~~~~~~~~~~~~~~ SLOT B ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_b <- subset(recordTableIndividual, DateTimeOriginal >= "2019-02-01 00:00:01" & DateTimeOriginal <= "2019-07-30 23:59:59")

## Generate capthist object while considering sampling effort 
CT_b <- subset(grid_traps, Station %in% slot_b$Station) # get CT table for slot B
camop_b <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_b))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_b))` >= "2019-02-01 00:00:01" & `as.Date(colnames(camop_b))` <= "2019-07-30 23:59:59")
dates <- colnames(camop_b) # get desired dates 
camop_b <- as.data.frame(t(camop_b)) # prepare matrix 
camop_b$date <- rownames(camop_b) # add dates as column for sub-setting
camop_b <- subset(camop_b, date %in% as.character(slot$`as.Date(colnames(camop_b))`)) # keep only dates in slot B
camop_b[,"date"] <- NULL # drop redundant dates
camop_b <- as.matrix(t(camop_b)) # rearrange as matrix
camop_b <- camop_b[rownames(camop_b) %in% slot_b$Station, ] # subset by stations found in record table 
slot_b$Species <- "jaguar"
# Build capthist object
cap_b <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_b,
                                           CTtable = CT_b,
                                           camOp = camop_b,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-02-01",
                                           includeEffort = TRUE)

## Examine capture history and trapping effort

summary(cap_b) # examine capthist
usage(traps(cap_b)) # examine trapping effort per occasion

plot(cap_b, tracks = TRUE, main = "Slot B") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_b))
hist(m, xlab = "Movement m", main = "Slot B")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_b, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_b, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_b), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first preliminary model
fit_b <- secr.fit(cap_b, buffer = 5*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_b # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_b, limits = FALSE) # plot detection probability

esa.plot(fit_b) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') # increase buffer slightly 
abline(v = 6 * initialsigma, lty = 2, col = 'green') # better

# Adjust mask buffer size
mask <- make.mask(traps(cap_b), buffer = 6*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_b.HN <- secr.fit(cap_b, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_b.EX <- secr.fit(cap_b, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_b.HN, EX = fit_b.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 6 * initialsigma)

# Final model slot B
fit_b <- secr.fit(cap_b, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")


## ~~~~~~~~~~~~~~~~~~~~~ SLOT C ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_c <- subset(recordTableIndividual, DateTimeOriginal >= "2019-03-01 00:00:01" & DateTimeOriginal <= "2019-08-27 23:59:59")

## Generate capthist object while considering sampling effort 
CT_c <- subset(grid_traps, Station %in% slot_c$Station) # get CT table for slot C
camop_c <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_c))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_c))` >= "2019-03-01 00:00:0" & `as.Date(colnames(camop_c))` <= "2019-08-27 23:59:59")
dates <- colnames(camop_c) # get desired dates 
camop_c <- as.data.frame(t(camop_c)) # prepare matrix 
camop_c$date <- rownames(camop_c) # add dates as column for sub-setting
camop_c <- subset(camop_c, date %in% as.character(slot$`as.Date(colnames(camop_c))`)) # keep only dates in slot C
camop_c[,"date"] <- NULL # drop redundant dates
camop_c <- as.matrix(t(camop_c)) # rearrange as matrix
camop_c <- camop_c[rownames(camop_c) %in% slot_c$Station, ] # subset by stations found in record table 
slot_c$Species <- "jaguar"
# Build capthist object
cap_c <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_c,
                                           CTtable = CT_c,
                                           camOp = camop_c,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-03-01",
                                           includeEffort = TRUE)

## Examine capture history and trapping effort

summary(cap_c) # examine capthist
usage(traps(cap_c)) # examine trapping effort per occasion

plot(cap_c, tracks = TRUE, main = "Slot C") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_c))
hist(m, xlab = "Movement m", main = "Slot C")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_c, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_c, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_c), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_c <- secr.fit(cap_c, buffer = 5*initialsigma, trace=FALSE, mask = mask, 
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_c # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_c, limits = FALSE) # plot detection probability

esa.plot(fit_c) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') # increase buffer slightly 
abline(v = 6 * initialsigma, lty = 2, col = 'green') # better

# Adjust mask buffer
mask <- make.mask(traps(cap_c), buffer = 6*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_c.HN <- secr.fit(cap_c, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_c.EX <- secr.fit(cap_c, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_c.HN, EX = fit_c.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 6 * initialsigma)

# Final model slot C
fit_c <- secr.fit(cap_c, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")

## ~~~~~~~~~~~~~~~~~~~~~ SLOT D ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_d <- subset(recordTableIndividual, DateTimeOriginal >= "2019-04-01 00:00:01" & DateTimeOriginal <= "2019-09-27 23:59:59")

## Generate capthist object while considering sampling effort 
CT_d <- subset(grid_traps, Station %in% slot_d$Station) # get CT table for slot D
camop_d <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_d))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_d))` >= "2019-04-01 00:00:01" & `as.Date(colnames(camop_d))` <= "2019-09-27 23:59:59")
dates <- colnames(camop_d) # get desired dates 
camop_d <- as.data.frame(t(camop_d)) # prepare matrix 
camop_d$date <- rownames(camop_d) # add dates as column for sub-setting
camop_d <- subset(camop_d, date %in% as.character(slot$`as.Date(colnames(camop_d))`)) # keep only dates in slot D
camop_d[,"date"] <- NULL # drop redundant dates
camop_d <- as.matrix(t(camop_d)) # rearrange as matrix
camop_d <- camop_d[rownames(camop_d) %in% slot_d$Station, ] # subset by stations found in record table 
slot_d$Species <- "jaguar"
# Build capthist object
cap_d <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_d,
                                           CTtable = CT_d,
                                           camOp = camop_d,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-04-01",
                                           includeEffort = TRUE)

## Examine capture history and trapping effort

summary(cap_d) # examine capthist
usage(traps(cap_d)) # examine trapping effort per occasion

plot(cap_d, tracks = TRUE, main = "Slot D") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_d))
hist(m, xlab = "Movement m", main = "Slot D")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_d, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_d, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_d), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_d <- secr.fit(cap_d, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_d # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_d, limits = FALSE) # plot detection probability

esa.plot(fit_d) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') 
abline(v = 6 * initialsigma, lty = 2, col = 'red') # better

# Adjust mask buffer size
mask <- make.mask(traps(cap_d), buffer = 6*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_d.HN <- secr.fit(cap_d, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_d.EX <- secr.fit(cap_d, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_d.HN, EX = fit_d.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 6 * initialsigma)

# Final model slot D
fit_d <- secr.fit(cap_d, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")

## ~~~~~~~~~~~~~~~~~~~~~ SLOT E ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_e <- subset(recordTableIndividual, DateTimeOriginal >= "2019-05-01 00:00:01" & DateTimeOriginal <= "2019-10-27 23:59:59")

## Generate capthist object while considering sampling effort 
CT_e <- subset(grid_traps, Station %in% slot_e$Station) # get CT table for slot E
camop_e <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_e))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_e))` >= "2019-05-01 00:00:01" & `as.Date(colnames(camop_e))` <= "2019-10-27 23:59:59")
dates <- colnames(camop_e) # get desired dates 
camop_e <- as.data.frame(t(camop_e)) # prepare matrix 
camop_e$date <- rownames(camop_e) # add dates as column for sub-setting
camop_e <- subset(camop_e, date %in% as.character(slot$`as.Date(colnames(camop_e))`)) # keep only dates in slot E
camop_e[,"date"] <- NULL # drop redundant dates
camop_e <- as.matrix(t(camop_e)) # rearrange as matrix
camop_e <- camop_e[rownames(camop_e) %in% slot_e$Station, ] # subset by stations found in record table 
slot_e$Species <- "jaguar"
# Build capthist object
cap_e <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_e,
                                           CTtable = CT_e,
                                           camOp = camop_e,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-05-01",
                                           includeEffort = TRUE)


## Examine capture history and trapping effort

summary(cap_e) # examine capthist
usage(traps(cap_e)) # examine trapping effort per occasion

plot(cap_e, tracks = TRUE, main = "Slot E") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_e))
hist(m, xlab = "Movement e", main = "Slot E")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_e, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_e, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_e), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_e <- secr.fit(cap_e, buffer = 5*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_e # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_e, limits = FALSE) # plot detection probability

esa.plot(fit_e) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') 
abline(v = 6 * initialsigma, lty = 2, col = 'red') # better

# Adjust mask buffer
mask <- make.mask(traps(cap_e), buffer = 6*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_e.HN <- secr.fit(cap_e, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_e.EX <- secr.fit(cap_e, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_e.HN, EX = fit_e.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 6 * initialsigma)

# Final model slot E
fit_e <- secr.fit(cap_e, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")

## ~~~~~~~~~~~~~~~~~~~~~ SLOT F ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_f <- subset(recordTableIndividual, DateTimeOriginal >= "2019-06-01 00:00:01" & DateTimeOriginal <= "2019-11-27 23:59:59")

## Generate capthist object while considering sampling effort 
CT_f <- subset(grid_traps, Station %in% slot_f$Station) # get CT table for slot \F
camop_f <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_f))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_f))` >= "2019-06-01 00:00:01" & `as.Date(colnames(camop_f))` <= "2019-11-28 23:59:59")
dates <- colnames(camop_f) # get desired dates 
camop_f <- as.data.frame(t(camop_f)) # prepare matrix 
camop_f$date <- rownames(camop_f) # add dates as column for sub-setting
camop_f <- subset(camop_f, date %in% as.character(slot$`as.Date(colnames(camop_f))`)) # keep only dates in slot F
camop_f[,"date"] <- NULL # drop redundant dates
camop_f <- as.matrix(t(camop_f)) # rearrange as matrix
camop_f <- camop_f[rownames(camop_f) %in% slot_f$Station, ] # subset by stations found in record table 
slot_f$Species <- "jaguar"
# Build capthist object
cap_f <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_f,
                                           CTtable = CT_f,
                                           camOp = camop_f,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-06-01",
                                           includeEffort = TRUE)


## Examine capture history and trapping effort

summary(cap_f) # examine capthist
usage(traps(cap_f)) # examine trapping effort per occasion

plot(cap_f, tracks = TRUE, main = "Slot F") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_f))
hist(m, xlab = "Movement f", main = "Slot F")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_f, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_f, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_f), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_f <- secr.fit(cap_f, buffer = 5*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_f # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_f, limits = FALSE) # plot detection probability

esa.plot(fit_f) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') 
abline(v = 8 * initialsigma, lty = 2, col = 'red') # better

# Adjust mask buffer
mask <- make.mask(traps(cap_f), buffer = 8*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_f.HN <- secr.fit(cap_f, buffer = 8*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_f.EX <- secr.fit(cap_f, buffer = 8*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_f.HN, EX = fit_f.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 10 * initialsigma)

# Final model slot F
fit_f <- secr.fit(cap_f, 
                  buffer = 8*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")

## ~~~~~~~~~~~~~~~~~~~~~ SLOT G ~~~~~~~~~~~~~~~~~~~~~

## Identify slot
slot_g <- subset(recordTableIndividual, DateTimeOriginal >= "2019-07-01 00:00:01" & DateTimeOriginal <= "2019-12-27 23:59:59")

## Generate capthist object while considering sampling effort 
CT_g <- subset(grid_traps, Station %in% slot_g$Station) # get CT table for slot \G
camop_g <- as.data.frame(camop_problem) # get camera operation table
slot <- as.data.frame(as.Date(colnames(camop_g))) # get slot
slot <- subset(slot, `as.Date(colnames(camop_g))` >= "2019-07-01 00:00:01" & `as.Date(colnames(camop_g))` <= "2019-12-27 23:59:59")
dates <- colnames(camop_g) # get desired dates 
camop_g <- as.data.frame(t(camop_g)) # prepare matrix 
camop_g$date <- rownames(camop_g) # add dates as column for sub-setting
camop_g <- subset(camop_g, date %in% as.character(slot$`as.Date(colnames(camop_g))`)) # keep only dates in slot G
camop_g[,"date"] <- NULL # drop redundant dates
camop_g <- as.matrix(t(camop_g)) # rearrange as matrix
camop_g <- camop_g[rownames(camop_g) %in% slot_g$Station, ] # subset by stations found in record table 
slot_g$Species <- "jaguar"
# Build capthist object
cap_g <- camtrapR::spatialDetectionHistory(recordTableIndividual = slot_g,
                                           CTtable = CT_g,
                                           camOp = camop_g,
                                           stationCol = "Station",
                                           speciesCol = "Species",
                                           species = "jaguar",
                                           output = "binary",
                                           recordDateTimeCol = "DateTimeOriginal",
                                           recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                           occasionLength = 1,
                                           individualCol = "ID",
                                           Xcol = "x",
                                           Ycol = "y",
                                           timeZone = "America/La_Paz",
                                           day1 = "2019-07-01",
                                           includeEffort = TRUE)


## Examine capture history and trapping effort

summary(cap_g) # examine capthist
usage(traps(cap_g)) # examine trapping effort per occasion

plot(cap_g, tracks = TRUE, main = "Slot G") # plot trapping history

## >>>>> Build SECR model <<<<< 

## Check movement distribution
m <- unlist(moves(cap_g))
hist(m, xlab = "Movement g", main = "Slot G")

# Check quick and biased population size, detection probability and home range estimates 
pfn(cap_g, N.estimator = c("n", "null","zippin","jackknife") )

# Estimate spatial scale parameter sigma
initialsigma <- RPSV(cap_g, CC = TRUE)
cat("Estimate of sigma =", initialsigma, "m\n")

# Build mask object
mask <- make.mask(traps(cap_g), buffer = 5*initialsigma, type = "trapbuffer")

# Fit a first, preliminary model
fit_g <- secr.fit(cap_g, buffer = 5*initialsigma, trace=FALSE, mask = mask,
                  start=list(D=0.0001, g0=0.1, sigma=1000))

fit_g # examine preliminary model 

## --- Check model diagnostics ---

plot(fit_g, limits = FALSE) # plot detection probability

esa.plot(fit_g) # check buffer size, it should reach a plateau 

abline(v = 5 * initialsigma, lty = 2, col = 'red') 
abline(v = 6 * initialsigma, lty = 2, col = 'red') # better

# Adjust mask buffer
mask <- make.mask(traps(cap_g), buffer = 6*initialsigma, type = "trapbuffer")

# ----- Choosing best detection function ------

# Half-normal
fit_g.HN <- secr.fit(cap_g, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "HN")
# Negative exponential
fit_g.EX <- secr.fit(cap_g, buffer = 6*initialsigma, trace=FALSE, mask = mask,
                     start=list(D=0.0001, g0=0.1, sigma=1000), detectfn = "EX")

fits <- secrlist(HN = fit_g.HN, EX = fit_g.EX) # list models together

predict(fits) # check estimates
AIC(fits)     # check AIC 

esa.plot(fits, max.buffer = 7 * initialsigma)

# Final model slot E
fit_g <- secr.fit(cap_g, 
                  buffer = 6*initialsigma, 
                  mask = mask,
                  trace=FALSE, 
                  start=list(D=0.0001, g0=0.1, sigma=1000),
                  detectfn = "HN")


## ------------------- 3.6 Get results & Draw Conclusions  --------------------

# List models
secr_grid <- secrlist(a = fit_a, b = fit_b, c = fit_c, d = fit_d, e = fit_e, f = fit_f, g = fit_g)

# Get density estimates
a <- predict(fit_a)[1, ]
b <- predict(fit_b)[1, ]
c <- predict(fit_c)[1, ]
d <- predict(fit_d)[1, ]
e <- predict(fit_e)[1, ]
f <- predict(fit_f)[1, ]
g <- predict(fit_g)[1, ]

## Combine estimates
secr_est <- as.data.frame(rbind(a,b,c,d,e,f,g))[ , -c(1,3)]
secr_est <- secr_est*10000 # transform to jaguars per 100 sq km 
rownames(secr_est) <- c("a", "b", "c", "d", "e", "f", "g")
secr_est$Session <- c(1:nrow(secr_est))

## Plot densities
dens_plot <- ggplot(data = secr_est) +
  geom_errorbar(aes(x = Session, ymin = lcl, ymax = ucl),
                alpha = .5, color = "grey40", width = .4, size = .3) +
  geom_line(aes(x = Session, y = estimate),
            linetype = "dashed", alpha = .5, color = "grey40", size = .3)+
  geom_point(aes(x = Session, y = estimate),
             size = 1.5, fill = "black", color = "black")+
  theme_bw(base_line_size = .2)+
  theme(axis.title.x=element_text(family="Arial", size = 6),
        axis.title.y=element_text(family="Arial", size = 6), 
        axis.text.x=element_text(family="Arial", size = 6),
        axis.text.y=element_text(family="Arial", size = 6))+
  ylab(bquote('Jaguars per 100 km'^2)) +
  xlab("Session")

# Export plot
ggsave(filename = "fig_4.png",
       plot = dens_plot,
       bg = "white",
       width = 80,
       height = 50,
       units = "mm",
       dpi = 1257)

## Get SECR parameters in a separated table
secr_final <- as.data.frame(cbind(
  rep(c("a", "b", "c", "d", "e", "f", "g"), each = 3),
  rep(c("D", "g0", "sigma"), times = 7)))
secr_final <- as.data.frame(cbind(secr_final, rbind(predict(fit_a),
                                                    predict(fit_b),
                                                    predict(fit_c),
                                                    predict(fit_d),
                                                    predict(fit_e),
                                                    predict(fit_f),
                                                    predict(fit_g))))
secr_final$link <- NULL
rownames(secr_final) <- NULL
colnames(secr_final) <- c("Session", "Value", "Estimate", "Std. Error", "Lower CI", "Upper CI")
secr_final[c(1,4,7,10,13,16,19), 3:6] <- secr_final[c(1,4,7,10,13,16,19), 3:6] * 10000 # get D per 100 sq km
secr_final[ ,c(3:6)] <- round(secr_final[ ,c(3:6)], digits = 3)

## Write final estimate table 
write.csv(secr_final, file = "secr_final.csv", row.names = FALSE)

# Clean global environment 
keep(act_plot, camtraps, dat.all, dat.clean, secr_grid, dens_plot, secr_final, recordTableIndividual, sure = TRUE)


################################# (4) Mapping ################################

## Scatter pie map using all data collected scaled by abundance relative to sampling effort

## Get temporally independent captures 
require(camtrapRdeluxe, quietly = TRUE, warn.conflicts = FALSE) 
recordTableIndividual <- assessTemporalIndependence(intable = dat.clean, 
                                                    deltaTimeComparedTo = "lastIndependentRecord", 
                                                    columnOfInterest    = "ID",
                                                    stationCol          = "Station", 
                                                    camerasIndependent  = FALSE, 
                                                    minDeltaTime        = 60)
detach("package:camtrapRdeluxe", unload = TRUE)

recordTableIndividual$Station <- str_trim(recordTableIndividual$Station, side = "both")
recordTableIndividual[is.na(recordTableIndividual$ID)]$ID <- "unidentified"

# Define factor levels of individual ID's
# recordTableIndividual$ID <- factor(recordTableIndividual$ID, 
#                                    levels = rev(c("Not identified",
#                                               "M-01","M-02","M-03","M-04","M-05","M-06",
#                                               "F-01","F-02","F-03","F-04",
#                                               "J-01-1","J-02-1","J-03-1","J-03-2","J-03-3")))
recordTableIndividual$ID <- factor(recordTableIndividual$ID, 
                                   levels = rev(c("M-01", "M-02", "F-03", "F-01", "J-03-1", "J-02-1", 
                                                  "J-03-2", "F-02", "M-05", "J-01-1", "M-04", "F-04",
                                                  "M-06", "J-03-3", "M-03", "unidentified")))

## Count unique individual captures per station 
individuals <- recordTableIndividual %>% group_by(ID) %>% 
  summarise(Count = n()) %>% arrange(ID)
count <- individuals$Count

# Set colour palette 
indPalette_c <- rev(c("#8b4513", "#008000", "#bdb76b", 
                  "#4b0082", "tomato3", "#48d1cc", "#ff4500", 
                  "#b03060", "#FCB714", "#00fa9a", "#0000ff", 
                  "#d8bfd8", "#ff00ff", "#1e90ff", "#ee82ee", "grey60"))
show_col(indPalette_c) # check colours 

# Plot Number of interdependent individual captures 
indi_bar <- ggplot(individuals, aes(x=ID, y = Count, fill=ID))+ # set input data 
  # set data format & visualization statistics 
  geom_bar(position = "dodge", stat = "identity", color = "transparent", alpha = 1)+
  # Set colour palette
  scale_fill_manual(values = indPalette_c) +
  geom_text(aes(label = Count), hjust = -.25)+
  # set themes
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(face = 1, size = 3, colour = "black",
                            family="sans"),
        legend.title = element_blank(),
        legend.text = element_text(face = 1, size = 8, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10),
                                   family="sans"),
        plot.title = element_text(hjust = 0.0, face = 2,family="sans", size = 6),
        axis.text.y = element_text(face = 1, size = 8, colour = "black",
                                   margin = margin(t = 0, r = 3, b = 0, l = 10),
                                   family="sans"),
        axis.text.x = element_text(face = 1, size = 7, colour = "black",
                                   margin = margin(t = 3, r = 0, b = 5, l = 0),
                                   family="sans"),
        axis.title.y = element_text(face = 1, size = 8, family="sans"),
        axis.title.x = element_text(face = 1, size = 8, family="sans"),
        axis.ticks = element_line(colour = "black"))+
  # set labels & extends 
  xlab("")+
  ggtitle(" ")+
  ylab("Captures")+
  expand_limits(y = 155)+
  guides(fill = guide_legend(reverse=T))+
  # flip x and y axes
  coord_flip()

## >>> Prepare data for the scatter pie map <<<

## Prepare input data
individuals <- dplyr::select(recordTableIndividual, lat, long, ID, Station)
ind.melt <- reshape2::melt(individuals, id.vars = c("lat", "long", "Station"))
ind.cast <- reshape2::dcast(ind.melt, Station~value, length) # ignore warning message
ind.cast$total <- rowSums(ind.cast[, -1])

## Get a DF with associated lat/long and Station values to add in
latlong <- dplyr::select(recordTableIndividual, lat, long, Station)
latlong <- unique(latlong)

## Merge data 
ind.pie <- merge(x=ind.cast,y=latlong,by.x="Station",by.y="Station")
dat.ind <- as.data.frame(table(recordTableIndividual$ID))
as.character(dat.ind$Var1)

## Number of individuals as label 
ind.pie$indis <- rowSums(ind.pie[ ,c(2:16)] >= 1)
ind.pie$total_div <- paste0(ind.pie$total, " (",ind.pie$indis, ")")

## Check if overlap is  a problem
plot(ind.pie$long, ind.pie$lat) # yes 
mapview(st_as_sf(ind.pie[ ,c("Station", "lat", "long")],
                 coords = c("long", "lat"), crs = 4326), legend = TRUE) 

# Move La Cachuela 
ind.pie[3, ]$long  <- ind.pie[3, ]$long + 0.002
ind.pie[3, ]$lat   <- ind.pie[3, ]$lat  - 0.002
# Move Orquidia 
ind.pie[6, ]$lat   <- ind.pie[6, ]$lat - 0.0035
ind.pie[6, ]$long  <- ind.pie[6, ]$long - 0.004
# Move G-02 
ind.pie[10, ]$long  <- ind.pie[10, ]$long + 0.0045
ind.pie[10, ]$lat   <- ind.pie[10, ]$lat - 0.004
# Move CaminoCachuela
ind.pie[1, ]$long <- ind.pie[1, ]$long + 0.0055
ind.pie[1, ]$lat  <- ind.pie[1, ]$lat  - 0.0013 
# Move G-01 
ind.pie[9, ]$lat   <- ind.pie[9, ]$lat + 0.0033
ind.pie[9, ]$long  <- ind.pie[9, ]$long - 0.002 
# Move Lagarto
ind.pie[4, ]$long  <- ind.pie[4, ]$long - 0.001 
ind.pie[4, ]$lat   <- ind.pie[4, ]$lat - 0.0015 
# Move Jaguar
ind.pie[2, ]$long  <- ind.pie[2, ]$long + 0.004
ind.pie[2, ]$lat   <- ind.pie[2, ]$lat + 0.0025
# Move Vaca Muerte
ind.pie[8, ]$long  <- ind.pie[8, ]$long + 0.0015 
# Move Laja de Noviquia
ind.pie[5, ]$lat  <- ind.pie[5, ]$lat - 0.002 
# Move G-05 North 
ind.pie[13, ]$lat   <- ind.pie[13, ]$lat + 0.002 

## Check again
plot(ind.pie$long, ind.pie$lat) # better
mapview(st_as_sf(ind.pie[ ,c("Station", "lat", "long")],
                 coords = c("long", "lat"), crs = 4326), legend = FALSE) 

## Get relative abundance index as scaling parameter 

## Calculate active days per Station
act_days <- as.data.frame(t(rbind(camtraps$Station, camtraps$Setup_date, camtraps$Retrieval_date))) # get setup and retrieval dates
act_days$V2 <- as.POSIXct(act_days$V2, format = "%d.%m.%Y") # set format 
act_days$V3 <- as.POSIXct(act_days$V3, format = "%d.%m.%Y") # set format 
act_days$difftime <- round(as.numeric(difftime(act_days$V3, act_days$V2, units = "days"))) # calculate days between setup and retrieval
act_days <- aggregate(difftime ~ V1, act_days, FUN = max, drop = FALSE) # keep only maximum per Station 
colnames(act_days) <- c("Station", "days_active")

## Count captures per Station 
station_count <- dplyr::select(recordTableIndividual %>% dplyr::group_by(Station) %>% 
                                 dplyr::summarise(Captures = n()), Captures, Station)

## Merge input data by station
rai <- na.omit(Hmisc::Merge(as.data.frame(station_count), as.data.frame(act_days), id = ~ Station))
row.names(rai) <- NULL

## calculate relative abundance index (see Bots et al.)
rai$index <- (rai$Captures / rai$days_active) * 1000
ind.pie$rai <- rai$index

## Get empty stations as an additional layer 
no_captures <- dplyr::select(subset(camtraps, !(Station %in% recordTableIndividual$Station)), Station, lat, long, NewGrid)
no_captures <- subset(no_captures, NewGrid == "y")
no_captures <- distinct(no_captures)
no_captures$lat <- as.numeric(no_captures$lat)
no_captures$long <- as.numeric(no_captures$long)

## ~~~~~~ Get base-map & Build scatter pie map ~~~~~~ 
register_google(key = "...") # Register Google API to access base maps 
basemap <- get_map(location = c(left = -62.0990, bottom = -16.3850, right = -61.8740, top = -16.3050), 
                   maptype = "satellite", source = "google", scale = 4) # Get base map
basemap_ref <- get_map(location = c(left = -62.0990, bottom = -16.3850, right = -61.8740, top = -16.3050), 
                   maptype = "satellite", source = "google", urlonly = TRUE) # get base map URL

## Categorize relative abundance index as a scaling parameter
ggplot(ind.pie, aes(x=rai)) + geom_histogram(binwidth = 20, color = "white")
ind.pie$rai <- as.numeric(ind.pie$rai) 

ind.pie$cat[ind.pie$rai <= 10] <- 0.0025                   
ind.pie$cat[ind.pie$rai > 10 & ind.pie$rai <= 20] <- 0.003
ind.pie$cat[ind.pie$rai > 20 & ind.pie$rai <= 40] <- 0.004 
ind.pie$cat[ind.pie$rai > 40 & ind.pie$rai <= 60] <- 0.005  
ind.pie$cat[ind.pie$rai > 60 & ind.pie$rai <= 80] <- 0.006  
ind.pie$cat[ind.pie$rai > 80 ] <- 0.007

## Build scatter pie map 
indi_pie <- 
  # set base-map & CRS projection
  ggmap(basemap) + coord_sf(crs = 4326) + 
  # set data input & format 
  geom_scatterpie(data = ind.pie, 
                  aes(long, lat, r = cat),
                  alpha = 1,
                  cols = as.character(dat.ind$Var1), size = .1,
                  color = "white") +
  # add empty stations 
  geom_spatial_point(data = no_captures, aes(long, lat), 
                     shape = 13, color = "grey96", size = 5,
                     inherit.aes = FALSE, crs = 4326) +
  # set color palette & thematic aesthetics 
  scale_fill_manual(values = indPalette_c)+
  # add empty stations 
  theme_linedraw()+ 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.0, face = 2, family="sans"),
        axis.text.y = element_text(face = 1, size = 6, family="sans", hjust = -1),
        axis.text.x = element_text(face = 1, size = 6, family="sans", vjust = -1),
        axis.title.y = element_text(face = 1, size = 8, family="sans"),
        axis.title.x = element_text(face = 1, size = 8, family="sans"),
        axis.ticks.length=unit(-0.15, "cm"),
        axis.ticks = element_line(colour = "white", size = .5))+
  # add north arrow
  # ggspatial::annotation_north_arrow(
  #   location = "tr", which_north = "true",
  #   style = ggspatial::north_arrow_nautical(
  #     fill = c("black", "white"),
  #     line_col = "black",
  #     text_family = "sans",
  #     text_col = "white")) +
  # add scale bar
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("black", "white"),
    text_family = "sans",
    pad_y = unit(.45, "cm"),
    text_col = "white") +
  # define labs 
  labs("") +
  xlab(" ")+
  ylab(" ") 

# Add bar plot as a legend & set title
jaguar_map <- ggarrange(indi_pie, indi_bar, ncol = 2, nrow = 1, widths = c(1.7,1.2))

# Export SECR map 
ggsave(filename = "Fig_5.png",
       plot = jaguar_map,
       dpi = 1000,
       height = 100,
       width = 170,
       units = "mm",
       bg = "white")

# Clean global environment
keep(act_plot, camtraps, dat.all, dat.clean, secr_final, jaguar_map, dens_plot, sure = TRUE)

################################# (5) Land cover analysis  ################################

# Read and prepare land cover data 
landuse <- read.csv("land_2019.csv", header = TRUE, sep = ",")
rownames(landuse) <-landuse$Type
landuse$Type <- NULL
landuse <- as.data.frame(t(landuse))
landuse$Station <- as.character(rownames(landuse))
landuse$Station <- gsub("\\.", "-", landuse$Station)
landuse[is.na(landuse)] <- 0
rownames(landuse) <- NULL

##  Prepare capture data
captures <- dat.all[, c(6,15,12,13,24,25,26,27,23,22)]   # relevant columns 
colnames(captures) <- c("ID", "DateTimeOriginal", "Age", "Sex", "Station", "Camera", "lat", "long","Grid", "Setup")
setDT(captures)
captures <- na.omit(captures, cols = c("DateTimeOriginal", "Station"), invert = FALSE)
DateTimeOriginal <- as.POSIXct(captures$DateTimeOriginal,
                               format = "%Y:%m:%d %H:%M:%S")
captures <- as.data.frame(apply(captures,2,function(x)gsub('\\s+', '',x)))
captures$DateTimeOriginal <- DateTimeOriginal 
captures <- subset(captures, Grid == "Grid" & 
                     DateTimeOriginal >= "2019-01-01 00:00:00" &
                     DateTimeOriginal <= "2019-12-31 23:59:59")
require(camtrapRdeluxe)
captures <- assessTemporalIndependence(intable = captures, # temp. independence  
                                       deltaTimeComparedTo = "lastIndependentRecord", 
                                       columnOfInterest    = "ID",
                                       stationCol          = "Station", 
                                       camerasIndependent  = FALSE, 
                                       minDeltaTime        = 60)
## Count captures per grid station
captures <- captures %>%  group_by(Station) %>% summarise(captures = n())
detach("package:camtrapRdeluxe", unload = TRUE)

##  Merge data  
captures$Station <- str_trim(captures$Station, side ="both")
landuse$Station <- str_trim(landuse$Station, side ="both")
landuse <- Hmisc::Merge(landuse, captures, all = TRUE, id = ~ Station)
landuse[is.na(landuse)] <- 0
colnames(landuse) <- c("station", "no_veg", "bare", "low_veg", "medium_veg", "high_veg", "captures")
# Add land use data based on cattle observations
landuse$landuse <- c("forest", "pasture", "pasture", "forest", "pasture", 
                     "forest", "forest", 
                     "pasture", "pasture", "pasture", "pasture")

## Test if captures differed significantly between forest and pasture
landuse %>% group_by(landuse) %>% get_summary_stats(captures, type = "median_iqr")

# Prepare boxplot
landuse_bxp <- ggplot(data = landuse, aes(x = landuse, y = captures)) +
  geom_boxplot(outlier.colour="transparent", outlier.size=0) +
  geom_jitter(shape=16, position=position_jitter(0.15)) +
  theme(axis.title.x=element_text(family="Arial", size = 6),
        axis.title.y=element_text(family="Arial", size = 6), 
        axis.text.x=element_text(family="Arial", size = 6),
        axis.text.y=element_text(family="Arial", size = 6))+
  ylab("Captures") +
  xlab(NULL)

# Test if captures differed between forest and pasture
landuse.test <- landuse %>% 
  wilcox_test(captures ~ landuse) %>%
  add_significance()
landuse.test

# Calculate the effect size and 
landuse %>% wilcox_effsize(captures ~ landuse)

# Visualize test results 
stat.test <- landuse.test %>% add_xy_position(x = "landuse")

landuse_bxp <- ggplot(data = landuse, aes(x = landuse, y = captures)) +
  # Boxplot
  geom_boxplot(outlier.colour="transparent", outlier.shape=8) +
  geom_jitter(shape=16, position=position_jitter(0.25), color = "grey30") +
  theme_bw(base_line_size = .2)+
  theme(axis.title.x=element_text(family="Arial", size = 6),
        axis.title.y=element_text(family="Arial", size = 6), 
        axis.text.x=element_text(family="Arial", size = 6),
        axis.text.y=element_text(family="Arial", size = 6),
        plot.title = element_text(family="Arial", size = 6))+
  ylab("Captures") +
  xlab(NULL)+ 
  # Test results 
  labs(title = get_test_label(landuse.test, detailed = TRUE, type = "text"))

# Export plot
ggsave(filename = "landuse.png",
       plot = landuse_bxp,
       bg = "white",
       width = 82,
       height = 72,
       units = "mm",
       dpi = 1257)

#####################################################################################
######################################## END ########################################
#####################################################################################

citation("camtrapR")   # Camera traps
citation("secr")       # SECR Analysis
citation("tidyverse")  # Graphics and data manipulation
citation("activity")   # Activity curves 
citation("data.table") # Data manipulation
citation("lubridate")  # Dates 
citation("ggmap")      # Base maps
citation("ggspatial")  # Mapping
citation("ggsn")       # Mapping
citation("scatterpie") # Mapping
citation("geosphere")  # Mapping
citation("scales")     # Mapping
citation("gridExtra")  # Mapping 
citation("GISTools")   # Mapping
citation("sp")         # Mapping
citation("sf")         # Mappingp
citation("pryr")       # Computations
citation("viridis")    # Aesthetics
citation("ggpubr")     # Plot arrangements 
citation("miscTools")  # misc
citation("extrafont")  # fonts
citation("cowplot")    # plot grid
citation("Hmisc")      # misc
