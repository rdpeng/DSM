library(splines)
f <- lm(Ozone ~ Solar.R + Wind + ns(Temp, 2), data = airquality)
m <- new("exposureModel", model = f, exposure = "Wind")


library(tsModel)
library(splines)
d0 <- readRDS("inst/ny.rds")
g <- glm(death ~ dow + ns(date, 19 * 6) + ns(tmpd, 6) + ns(rmtmpd, 6) + ns(dptp, 3) + ns(rmdptp, 3) + Lag(pm10tmean, 1), data = d0, family = poisson)

m <- new("exposureModel", model = g, exposure = "pm10tmean")
ap <- newAPTS(g, "pm10tmean")

