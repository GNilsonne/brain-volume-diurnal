# Script for preregistration: diurnal variation of brain gray and white matter volume
# Gustav Nilsonne

# This script analyses metadata from datasets to be included, 
# builds models, visualises predictions, and performs simulations

# Require packages
require(nlme)
require(effects)

# Read files
# These files are metadata files from NITRC
# Some of the files were in .csv format and some were in .txt format

myFiles <- list.files(pattern = "*.csv")
myFiles <- myFiles[myFiles %in% c("maclaren_scan_number_vs_time_of_day.csv", "anatMRIQC_ds000030_R1.4.csv") == F] # Remove Maclaren data
for(i in 1:length(myFiles)){
  data <- read.csv(myFiles[i])
  data$filename <- myFiles[i]
  if(!exists("data_csvfiles")){
    data_csvfiles <- data[, c("SUBID", "SESSION", "SEX", "TIME_OF_DAY", "filename")]
  } else {
    data_csvfiles <- rbind(data_csvfiles, data[, c("SUBID", "SESSION", "SEX", "TIME_OF_DAY", "filename")])
  }
}

data_csvfiles$TIME_OF_DAY <- substr(data_csvfiles$TIME_OF_DAY, 1, 1) # Remove additional characters in MRN file

myFiles2 <- list.files(pattern = "*.txt")
myFiles2 <- myFiles2[myFiles2 != "russ_reduced.txt"] # Remove Poldrack data
for(i in 1:length(myFiles2)){
  data <- read.csv(myFiles2[i])
  data$filename <- myFiles2[i]
  if(!exists("data_txtfiles")){
    data_txtfiles <- data[, c("SUBID", "SESSION", "SEX", "TIME_OF_DAY", "filename")]
  } else {
    data_txtfiles <- rbind(data_txtfiles, data[, c("SUBID", "SESSION", "SEX", "TIME_OF_DAY", "filename")])
  }
}

data <- rbind(data_csvfiles, data_txtfiles)

# Find and include only participants that have scans at more than one time point
SUBID_unbalanced <- vector()
for(i in unique(data$SUBID)){
  tempdata <- data[data$SUBID == i, ]
  if(length(unique(tempdata$TIME_OF_DAY)) > 1 ){
    SUBID_unbalanced <- c(SUBID_unbalanced, i)
  }
}
data_unbalanced <- data[data$SUBID %in% SUBID_unbalanced, ]

data_unbalanced$dataset <- substr(data_unbalanced$filename,1,nchar(data_unbalanced$filename)-4)

# Tabulate and plot number of scans by time of day
counts <- table(data_unbalanced$dataset, data_unbalanced$TIME_OF_DAY)
counts
tol18rainbow = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
barplot(counts, names.arg = c("4-8", "8-12", "12-16", "16-20", "20-24"), xlab = "Time of day (bin)", ylab = "n scans", col = tol18rainbow, main = "Unique scans")
legend("topright", legend = rownames(counts), bty = "n", fill = tol18rainbow, cex = 0.8)

# Tabulate and plot number of scans by time of day for unique participants
no_scans <- colSums(table(data_unbalanced$SUBID, data_unbalanced$dataset))
no_scans
sum(no_scans)
data_for_number_of_participants <- data_unbalanced
duplicated_indices2 <- duplicated(data_unbalanced[, c("SUBID")])
data_unbalanced_uniqueSUBID2 <- data_unbalanced[duplicated_indices2 == FALSE, ]
number_of_participants <- colSums(table(data_unbalanced_uniqueSUBID2$SUBID, data_unbalanced_uniqueSUBID2$dataset))
number_of_participants
sum(number_of_participants)
counts2 <- table(data_unbalanced_uniqueSUBID2$dataset, data_unbalanced_uniqueSUBID2$TIME_OF_DAY)
barplot(counts2, names.arg = c("4-8", "8-12", "12-16", "16-20", "20-24"), xlab = "Time of day (bin)", ylab = "n participants", col = tol18rainbow, main = "Unique participants")
legend("topright", legend = rownames(counts), bty = "n", fill = tol18rainbow, cex = 0.8)

# Read additional data from Maclaren et al.
maclaren <- read.csv("maclaren_scan_number_vs_time_of_day.csv")
names(maclaren) <- c("scan_no", "s1", "s2", "s3")

maclaren$s1 <- as.character(maclaren$s1)
maclaren$s1[nchar(maclaren$s1) == 3] <- paste("0", maclaren$s1[nchar(maclaren$s1) == 3], sep = "")
maclaren$s1_hm <- as.POSIXct(x = as.character(maclaren$s1), format = "%H%M")
maclaren$s2 <- as.character(maclaren$s2)
maclaren$s2[nchar(maclaren$s2) == 3] <- paste("0", maclaren$s2[nchar(maclaren$s2) == 3], sep = "")
maclaren$s2_hm <- as.POSIXct(x = as.character(maclaren$s2), format = "%H%M")
maclaren$s3 <- as.character(maclaren$s3)
maclaren$s3[nchar(maclaren$s3) == 3] <- paste("0", maclaren$s3[nchar(maclaren$s3) == 3], sep = "")
maclaren$s3_hm <- as.POSIXct(x = as.character(maclaren$s3), format = "%H%M")

maclaren$s1_minutesfromseven <- as.numeric(substr(as.character(maclaren$s1), 1, 2))*60 + as.numeric(substr(as.character(maclaren$s1), 3, 4)) - 7*60
maclaren$s2_minutesfromseven <- as.numeric(substr(as.character(maclaren$s2), 1, 2))*60 + as.numeric(substr(as.character(maclaren$s2), 3, 4)) - 7*60
maclaren$s3_minutesfromseven <- as.numeric(substr(as.character(maclaren$s3), 1, 2))*60 + as.numeric(substr(as.character(maclaren$s3), 3, 4)) - 7*60

min <- min(c(maclaren$s1_hm, maclaren$s2_hm, maclaren$s3_hm))
max <- max(c(maclaren$s1_hm, maclaren$s2_hm, maclaren$s3_hm))

# Plot Maclaren dataset
plot(rep(1, 40) ~ maclaren$s1_hm, frame.plot = F, xlim = c(min, max), yaxt = "n", ylim = c(0, 4), xlab = "time", ylab = "", type = "n", main = "Maclaren", xaxt = "n")
abline(h = c(1, 2, 3), lty = 2, col = "gray")
points(rep(1, 40) ~ maclaren$s1_hm, pch = "|")
points(rep(2, 40) ~ maclaren$s2_hm, pch = "|")
points(rep(3, 40) ~ maclaren$s3_hm, pch = "|")
axis(1, at = as.POSIXct(x = c("08:00", "12:00", "16:00", "20:00"), format="%H:%M"), labels = c("08:00", "12:00", "16:00", "20:00"))
axis(2, type = "n", at = c(1, 2, 3), labels = c("S1", "S2", "S3"), tick = F, las = 1)

maclaren2 <- data.frame(SUBID = c(rep(1, 40), rep(2, 40), rep(3, 40)), SESSION = 1:40, SEX = NA, TIME_OF_DAY = NA, filename = NA, dataset = "Maclaren", tod = c(maclaren$s1_hm, maclaren$s2_hm, maclaren$s3_hm))

# Format time of day and merge datasets
data_unbalanced$tod <- as.POSIXct(x = "0700", format = "%H%M")
data_unbalanced$tod[data_unbalanced$TIME_OF_DAY == 2] <- as.POSIXct(x = "1000", format = "%H%M")
data_unbalanced$tod[data_unbalanced$TIME_OF_DAY == 3] <- as.POSIXct(x = "1400", format = "%H%M")
data_unbalanced$tod[data_unbalanced$TIME_OF_DAY == 4] <- as.POSIXct(x = "1800", format = "%H%M")
data_unbalanced$tod[data_unbalanced$TIME_OF_DAY == 5] <- as.POSIXct(x = "2100", format = "%H%M")

maclaren3 <- data.frame(SUBID = c(rep(1, 40), rep(2, 40), rep(3, 40)), SESSION = 1:40, SEX = NA, TIME_OF_DAY = NA, filename = NA, dataset = "Maclaren", tod = c(maclaren$s1_hm, maclaren$s2_hm, maclaren$s3_hm))

data_unbalanced <- rbind(data_unbalanced, maclaren2)

# Make new time format in minutes starting from 07:00
data_unbalanced$minutesfromseven <- -60 # Failsafe placeholder datapoint
data_unbalanced$minutesfromseven[data_unbalanced$TIME_OF_DAY == 1] <- 0
data_unbalanced$minutesfromseven[data_unbalanced$TIME_OF_DAY == 2] <- 180
data_unbalanced$minutesfromseven[data_unbalanced$TIME_OF_DAY == 3] <- 420
data_unbalanced$minutesfromseven[data_unbalanced$TIME_OF_DAY == 4] <- 660
data_unbalanced$minutesfromseven[data_unbalanced$TIME_OF_DAY == 5] <- 900

data_unbalanced$minutesfromseven[data_unbalanced$SUBID == 1] <- maclaren$s1_minutesfromseven
data_unbalanced$minutesfromseven[data_unbalanced$SUBID == 2] <- maclaren$s2_minutesfromseven
data_unbalanced$minutesfromseven[data_unbalanced$SUBID == 3] <- maclaren$s3_minutesfromseven

# Build model using simulated data
set.seed(1)
ds000030 <- read.csv("anatMRIQC_ds000030_R1.4.csv") # Example data from openfmri dataset ds000030 to provide resonable values
data_unbalanced$gm_sim <- rnorm(mean = mean(ds000030$icvs_gm), sd = sd(ds000030$icvs_gm), n = length(data_unbalanced$SUBID))
lm_sim_lin <- lme(gm_sim ~ minutesfromseven, data = data_unbalanced, random = ~ minutesfromseven|dataset/SUBID, control = lmeControl(opt = 'optim'))
summary(lm_sim_lin)
plot(lm_sim_lin)
plot(effect("minutesfromseven", lm_sim_lin), main = "Predicted, linear model", xlab = "minutes from 07:00", ylab = "gray matter fraction")
lm_sim_exp <- lme(gm_sim ~ log(1 + minutesfromseven), data = data_unbalanced, random = ~ log(1 + minutesfromseven)|dataset/SUBID, control = lmeControl(opt = 'optim'))
summary(lm_sim_exp)
plot(lm_sim_exp)
plot(effect("log(1 + minutesfromseven)", lm_sim_exp), main = "Predicted, exponential model", xlab = "minutes from 07:00", ylab = "gray matter fraction")
rm(.Random.seed, envir=globalenv()) # Reset random seed

x_min <- min(data_unbalanced$tod)
x_max <- max(data_unbalanced$tod)
y_min <- min(data_unbalanced$gm_sim)
y_max <- max(data_unbalanced$gm_sim)
for(i in 1:length(unique(data_unbalanced$dataset))){
  plot(gm_sim ~ tod, data = data_unbalanced[data_unbalanced$dataset == unique(data_unbalanced$dataset)[i], ], frame.plot = F, xlab = "time of day", ylab = "gray matter fraction", main = unique(data_unbalanced$dataset)[i], xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  abline(a = lm_sim$coefficients$fixed[1] + lm_sim$coefficients$random$dataset[i, 1], b = lm_sim$coefficients$fixed[2] + lm_sim$coefficients$random$dataset[i, 2])
}


# Test models on simulated data
n_iterations <- 10
sim_results <- data.frame(matrix(nrow = n_iterations, ncol = 15)) # Initialise data frame
names(sim_results) <- c("iteration", "lower_lin", "est_lin", "upper_lin", "p_lin", 
                        "logLik_lin", "AIC_lin", "lower_exp", "est_exp", "upper_exp", 
                        "p_exp", "logLik_exp", "AIC_exp", "p_exp", "p_lin_vs_exp")
pb <- winProgressBar(title = "progress bar", min = 0, max = n_iterations, width = 300)
for(i in 1:n_iterations){
  data_unbalanced$gm_sim <- rnorm(mean = mean(ds000030$icvs_gm), sd = sd(ds000030$icvs_gm), n = length(data_unbalanced$SUBID))
  lm_lin <- tryCatch(lme(gm_sim ~ minutesfromseven, data = data_unbalanced, random = ~ minutesfromseven|dataset/SUBID, control = lmeControl(opt = 'optim')), error = function(e) NA)
  lm_exp <- tryCatch(lme(gm_sim ~ log(1 + minutesfromseven), data = data_unbalanced, random = ~ log(1 + minutesfromseven)|dataset/SUBID, control = lmeControl(opt = 'optim')), error = function(e) NA)
  if(length(lm_lin) == 18 & length(lm_exp) == 18){
    out <- data.frame(iteration = i, 
                      lower_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 1],
                      est_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 2],
                      upper_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 3],
                      p_lin = summary(lm_lin)$tTable[2, 5], 
                      logLik_lin = logLik(lm_lin, REML = F),
                      AIC_lin = AIC(lm_lin),
                      lower_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 1],
                      est_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 2],
                      upper_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 3],
                      p_exp = summary(lm_exp)$tTable[2, 5], 
                      logLik_exp = logLik(lm_exp, REML = F),
                      AIC_exp = AIC(lm_exp),
                      p_exp = summary(lm_exp)$tTable[2, 5],
                      p_lin_vs_exp = 1 - pchisq(abs(logLik(lm_lin, REML = F)[1] - logLik(lm_exp, REML = F)[1])*2, 2))
  } else if(length(lm_lin) == 18){
    out <- data.frame(iteration = i, 
                      lower_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 1],
                      est_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 2],
                      upper_lin <- intervals(lm_lin, which = "fixed")$fixed[2, 3],
                      p_lin = summary(lm_lin)$tTable[2, 5], 
                      logLik_lin = logLik(lm_lin, REML = F),
                      AIC_lin = AIC(lm_lin),
                      lower_exp <- NA,
                      est_exp <- NA,
                      upper_exp <- NA,
                      p_exp = NA, 
                      logLik_exp = NA,
                      AIC_exp = NA,
                      p_exp = NA,
                      p_lin_vs_exp = NA)
  } else {
    out <- data.frame(iteration = i, 
                      lower_lin <- NA,
                      est_lin <- NA,
                      upper_lin <- NA,
                      p_lin = NA, 
                      logLik_lin = NA,
                      AIC_lin = NA,
                      lower_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 1],
                      est_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 2],
                      upper_exp <- intervals(lm_exp, which = "fixed")$fixed[2, 3],
                      p_exp = summary(lm_exp)$tTable[2, 5], 
                      logLik_exp = logLik(lm_exp, REML = F),
                      AIC_exp = AIC(lm_exp),
                      p_exp = summary(lm_exp)$tTable[2, 5],
                      p_lin_vs_exp = 1 - pchisq(abs(logLik(lm_lin, REML = F)[1] - logLik(lm_exp, REML = F)[1])*2, 2))
  }
  sim_results[i, ] <- out
  setWinProgressBar(pb, i, title=paste( round(i/n_iterations*100, 0), "% done"))
}
close(pb)

hist(sim_results$p_lin, breaks = 19, main = "Linear model: p-values for time of day")
hist(sim_results$p_exp, breaks = 19, main = "Exponential model: p-values for time of day")

d <- density(sim_results$logLik_lin)
plot(d, main = "Log likelihood", frame.plot = F, type = "n")
polygon(d, col="gray", border="gray")
lines(density(sim_results$logLik_exp), lwd = 2)

hist(sim_results$p_lin_vs_exp, breaks = 19, main = "Model comparison: p-valused for LogLik-test")
