################################################################################
## Data structures for environmental epi example
################################################################################

setGeneric("units")
setGeneric("units<-")

setOldClass("lm")
setOldClass("glm")

setClass("exposureModel",
         representation(model = "lm",
                        exposure = "character"))

setGeneric("getNames", function(x, ...) standardGeneric("getNames"))
setMethod("getNames", "exposureModel",
          function(x, ...) {
                  frm <- formula(x@model)
                  others <- attr(terms(x@model), "term.labels")
                  i <- grep(x@exposure, others)
                  list(response = as.character(frm)[2],
                       exposure = others[i],
                       other = others[-i])
          })

char2comma <- function(x) {
        n <- length(x)
        if(n <= 1L)
                x
        else if(n == 2L)
                paste(sQuote(x), collapse = " and ")
        else
                paste0(paste(sQuote(x[-n]), collapse = ", "), ", and ",
                       sQuote(x[n]))
}

setGeneric("modelType", function(object, ...) standardGeneric("modelType"))

setMethod("modelType", "lm",
          function(object, ...) {
                  "linear regression"
          })

setMethod("modelType", "glm",
          function(object, ...) {
                  fam <- family(object)
                  if(fam$family == "gaussian")
                          "linear regression"
                  else if(fam$family == "binomial") {
                          paste(fam$link, "regression")
                  }
                  else 
                          paste(fam$family, "regression")                  
          })

setGeneric("describeModel", function(object, ...)
           standardGeneric("describeModel"))

setMethod("describeModel", "exposureModel", function(object, ...) {
        nms <- getNames(object)
        sprintf("a %s model where %s was regressed on exposure %s along with %s",
                modelType(object@model), sQuote(nms$response),
                sQuote(nms$exposure), char2comma(nms$other))
})

setGeneric("modelEffect", function(object, ...) standardGeneric("modelEffect"))


modelEffect.method <- function(object, incr = 1, pvalue = FALSE, expcoef = TRUE, ...) {
        fam <- family(object@model)$family
        nms <- getNames(object)
        ci <- suppressMessages(confint(object@model, nms$exposure))
        ci <- ci * incr
        cf <- coef(object@model)
        i <- grep(object@exposure, names(cf))
        b <- cf[i]
        bb <- b * incr

        if(fam == "binomial") {
                bb <- exp(bb)
                ci <- exp(ci)
                ci.txt <- prettyNum(ci)
                msg <- sprintf("a %s unit %s in %s was associated with a %s (95%% CI: %s, %s) %s odds of %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), prettyNum(abs(bb)),
                               ci.txt[1], ci.txt[2],
                               ifelse(bb > 0, "increased", "decreased"),
                               sQuote(nms$response))
        }
        else if(fam == "poisson") {
                bb <- exp(bb)
                ci <- exp(ci)
                ci.txt <- prettyNum(ci)
                msg <- sprintf("a %s unit %s in %s was associated with a %s (95%% CI: %s, %s) %s rate of %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), prettyNum(abs(bb)),
                               ci.txt[1], ci.txt[2],
                               ifelse(bb > 0, "increased", "decreased"),
                               sQuote(nms$response))
        }
        else {
                ci.txt <- prettyNum(ci)
                msg <- sprintf("a %s unit %s in %s was associated with a %s (95%% CI: %s, %s) unit %s in %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), prettyNum(abs(bb)),
                               ci.txt[1], ci.txt[2],
                               ifelse(bb > 0, "increase", "decrease"),
                               sQuote(nms$response))
        }
        if(pvalue) {
                p <- summary(object@model)$coefficients[nms$exposure, 4]
                msg <- sprintf("%s (p = %s)", msg, prettyNum(p))
        }
        msg
}


setMethod("modelEffect", "exposureModel", modelEffect.method)

startcaps <- function(x) {
}

setMethod("show", "exposureModel",
          function(object) {
                  out <- describeModel(object)
                  cat(out, ".\n", sep = "")
                  invisible(object)
          })


setClass("airpollutionTS",
         representation(timevar = "character",
                        "exposureModel"))
newAPTS <- function(model, exposure, timevar = NULL) {
        expm <- new("exposureModel", model = model, exposure = exposure)
        if(is.null(timevar)) {
                nms <- getNames(expm)
                i <- grep(c("date|time"), tolower(nms$other), perl = TRUE)
                if(length(i) > 0) {
                        timevar <- nms$other[i]
                        message(sprintf("using '%s' for 'timevar'", timevar))
                }
                else
                        stop("need to specify 'timevar'")
        }
        new("airpollutionTS", timevar = timevar, expm)
}

setGeneric("sensitivityAnalysis", function(object, ...)
           standardGeneric("sensitivityAnalysis"))

setMethod("sensitivityAnalysis", "airpollutionTS", function(object, ...) {
        
})



adjustTimeDF <- function(object, dfseq, timeVar = "time", smoothType = "ns") {
        confounderFormula <- as.formula(object$call$confounders)

        origTimeVar <- grep(timeVar, attr(terms(confounderFormula), "term.labels"),
                            value = TRUE)
        stopifnot(length(origTimeVar) == 1)
        
        timeVec <- paste("ns(", timeVar, ", ", dfseq, ")", sep = "")
        results <- vector("list", length = length(dfseq))

        for(i in seq(along = dfseq)) {
                newTimeVar <- timeVec[i]
                newFormula <- as.formula(paste("~ . -", origTimeVar, "+", newTimeVar))
                results[[i]] <- update(object, confounders = newFormula)
        }
        if(length(results) == 1)
                results <- results[[1]]
        structure(results, class = "adjustTimeDF", dfseq = dfseq)
}
