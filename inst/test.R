library(splines)
f <- lm(Ozone ~ Solar.R + Wind + ns(Temp, 2), data = airquality)
m <- new("exposureModel", model = f, exposure = "Wind")


library(tsModel)
library(splines)
source("R/datastructures.R")
## d0 <- readAPTSData("inst/ny.rds")
d0 <- readAPTSData("inst/chic.rds")
ch <- checkData(d0)
g <- glm(death ~ dow + ns(date, 19 * 6) + ns(tmpd, 6) + ns(rmtmpd, 6) + ns(dptp, 3) + ns(rmdptp, 3) + Lag(pm10tmean, 1), data = as(d0, "data.frame"), family = poisson)
ap <- new("APTSModel", timevar = "date", model = g, exposure = "pm10tmean", summary = ch)

sens <- sensitivityAnalysis(ap, 19 * 2:12)

