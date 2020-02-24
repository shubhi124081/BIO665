source('/Users/zeyi/Documents/PhD/Classes/Bayesian_inference/JimClarkFunctions2020.R')

summ <- data.frame(mean = colMeans(chains), 
                   median = apply(chains,2,median), 
                   Quantile1 = apply(chains,2,quantile,c(.025,.975))[1,],
                   Quantile2 = apply(chains,2,quantile,c(.025,.975))[2,],
                   cor_Mean = cor(chains)[1,],
                   cor_Variance = cor(chains)[2,])

row.names(summ)[1] <- c("Mean")
row.names(summ)[2] <- c("Variance")


knitr::kable(
  summ, 
  caption = "Summary for Chains of Mean and Variance Parameter"
)
