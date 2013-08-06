# Benchmark analysis evaluating method for choosing the degrees of
# freedom for the smooth function of time


```r
library(splines)
library(tsModel)
```

```
## Time Series Modeling for Air Pollution and Health (0.6)
```


This function chooses the degrees of freedom for a given dataset by
predicting PM10. The output is an integer degrees of freedom per year.


```r
getDF <- function(dataset, dfseq = 2:20) {
    aic <- sapply(dfseq, function(dfuse) {
        dfy <- dfuse * 8L
        form <- reformulate(c("ns(tmpd, 6)", sprintf("ns(time, %d)", dfy)), 
            response = "pm10tmean")
        model <- lm(form, data = dataset)
        AIC(model)
    })
    dfseq[which.min(aic)]
}
```





Here are the final summarized results with bias, standard error, and
root mean squared error.


```r
print(final)
```

```
##                             Bias        SE      RMSE
## simDatasets-fg.rda   0.000011043 0.0002547 0.0002549
## simDatasets-fg10.rda 0.000009232 0.0008059 0.0008060
## simDatasets-gf.rda   0.000002076 0.0002587 0.0002587
## simDatasets-gf10.rda 0.000032615 0.0008182 0.0008188
```




