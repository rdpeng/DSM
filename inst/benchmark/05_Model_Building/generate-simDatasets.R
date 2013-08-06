######################################################################
## Functions for generating data

simFramework <- function(data, simControl, method) {
    df.true <- simControl$df.true
    df.truex <- simControl$df.truex
    if(is.null(simControl$agecat)) {
        ## Collapse age groups
        death <- rowSums(sapply(split(data, data$agecat), "[[", "death"))
    }
    else if(simControl$agecat == "over75") 
       death <- subset(data, agecat == 3)$death
    else if(simControl$agecat == "under75")
       death <- subset(data, agecat != 3)$death
    else
        stop("Incorrect age category specified")
    sdata <- data[data$agecat == 1, ]
    year <- trunc(sdata[, "date"] / 10000)
    use <- year <= simControl$maxYear
    time <- sdata[use, "time"]
    death <- death[use]
    tmpd <- sdata[use, "tmpd"]
    pm10tmean <- sdata[use, "pm10tmean"]
    pm10tmean[is.na(pm10tmean)] <- mean(pm10tmean, na.rm = TRUE)

    if(method == "glm") {
        library(splines)    
        mortFit <- glm(death ~ ns(time, df.true[1]) + ns(tmpd, df.true[2]),
                       family = poisson, na.action = na.exclude)
        pollFit <- lm(pm10tmean ~ ns(time, df.truex[1])
                      + ns(tmpd, df.truex[2]), na.action = na.exclude)
        sigma2.xi <- (summary(pollFit)$sigma)^2
    }
    else if(method == "gamS") {
        ## This must be running on S-PLUS!
        mortFit <- gam(death ~ s(time, df.true[1]) + s(tmpd, df.true[2]),
                       family = poisson, na.action = na.exclude,
                       control = gam.control(epsilon = 1e-10, maxit = 200,
                       bf.maxit = 200))
        pollFit <- gam(pm10tmean ~ s(time, df.truex[1])
                       + s(tmpd, df.truex[2]), na.action = na.exclude,
                       control = gam.control(epsilon = 1e-10, maxit = 200,
                       bf.maxit = 200))
        r <- residuals(pollFit)
        sigma2.xi <- sum(r^2)/(length(r) - sum(df.truex))
    }
    else if (method == "gamR") {
        ## This should be running in mgcv 1.0-5 or higher
        stop("Not finished yet")

        G <- gam(death ~ s(time, k = 50 * 8) + s(tmpd, k = 10),
                 sp = 1, fit = FALSE)
    }
    else
        stop("Wrong ", sQuote("method"), " supplied")
    
    list(mortFit = mortFit, pollFit = pollFit, tmpd = tmpd, time = time,
         sigma2.xi = sigma2.xi)
}

## Generate a simulated pollution/mortality dataset

simData <- function(simFramework, simControl) {
    varScale <- simControl$varScale
    trueBeta <- simControl$trueBeta
    if(varScale <= 0)
        stop("varScale <= 0")
    pMort <- predict(simFramework$mortFit)
    pPoll <- predict(simFramework$pollFit)

    if(length(pMort) != length(pPoll))
        stop("length(pMort) != length(pPoll)")
    if(any(is.na(pMort)))
        stop("any(is.na(pMort))")
    if(any(is.na(pPoll)))
        stop("any(is.na(pPoll))")

    n <- length(pMort)

    if(!is.null(simControl$noSimPoll))
        simPoll <- model.frame(simFramework$pollFit)$pm10tmean
    else
        simPoll <- rnorm(n, pPoll, sqrt(simFramework$sigma2.xi / varScale))
    simMort <- rpois(n, exp(trueBeta * simPoll + pMort))

    simdf <- data.frame(death = simMort, pm10tmean = simPoll,
                        tmpd = simFramework$tmpd, time = simFramework$time)
    attr(simdf, "beta") <- trueBeta
    simdf
}


genDatasets <- function(Nsim, origData, simControl,
                        method = c("glm", "gamS", "gamR")) {
    method <- match.arg(method)
    simF <- simFramework(origData, simControl, method)
    datasets <- vector("list", length = Nsim)

    for(i in seq(along = datasets)) {
        cat(i, "\n")
        datasets[[i]] <- simData(simF, simControl)
    }
    attr(datasets, "simControl") <- simControl
    datasets
}


######################################################################
## Generate 500 simulated datasets

library(NMMAPSdata)
loadCity("minn")

Nsim <- 500
maxYear <- 1994
nYears <- maxYear - 1987 + 1

## f more wiggly than g, low concurvity

simControl <- list(
                   df.true = c(7 * nYears, 6),
                   df.truex = c(4 * nYears, 3),
                   varScale = 1,  
                   trueBeta = 0,
                   maxYear = maxYear,
                   )
seed <- 1000
set.seed(seed)
simDatasets <- genDatasets(Nsim, minn, simControl)

attr(simDatasets, "simControl") <- simControl
attr(simDatasets, "seed") <- seed

dump("simDatasets", file = "simulated-datasets/simDatasets-fg.R")
save(simDatasets, file = "simulated-datasets/simDatasets-fg.rda",
     compress = TRUE)


## f more wiggly than g, high concurvity

simControl <- list(
                   df.true = c(7 * nYears, 6),
                   df.truex = c(4 * nYears, 3),
                   varScale = 10,  
                   trueBeta = 0,
                   maxYear = maxYear,
                   )
seed <- 2000
set.seed(seed)
simDatasets <- genDatasets(Nsim, minn, simControl)

attr(simDatasets, "simControl") <- simControl
attr(simDatasets, "seed") <- seed

dump("simDatasets", file = "simulated-datasets/simDatasets-fg10.R")
save(simDatasets, file = "simulated-datasets/simDatasets-fg10.rda",
     compress = TRUE)


## g more wiggly than f, low concurvity

simControl <- list(
                   df.true = c(4 * nYears, 3),
                   df.truex = c(7 * nYears, 6),
                   varScale = 1,  
                   trueBeta = 0,
                   maxYear = maxYear,
                   )
seed <- 3000
set.seed(seed)
simDatasets <- genDatasets(Nsim, minn, simControl)

attr(simDatasets, "simControl") <- simControl
attr(simDatasets, "seed") <- seed

dump("simDatasets", file = "simulated-datasets/simDatasets-gf.R")
save(simDatasets, file = "simulated-datasets/simDatasets-gf.rda",
     compress = TRUE)


## g more wiggly than f, low concurvity

simControl <- list(
                   df.true = c(4 * nYears, 3),
                   df.truex = c(7 * nYears, 6),
                   varScale = 10,  
                   trueBeta = 0,
                   maxYear = maxYear,
                   )
seed <- 4000
set.seed(seed)
simDatasets <- genDatasets(Nsim, minn, simControl)

attr(simDatasets, "simControl") <- simControl
attr(simDatasets, "seed") <- seed

dump("simDatasets", file = "simulated-datasets/simDatasets-gf10.R")
save(simDatasets, file = "simulated-datasets/simDatasets-gf10.rda",
     compress = TRUE)
