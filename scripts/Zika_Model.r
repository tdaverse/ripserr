devtools::load_all()
data(aegypti)
data(zikaPredict)

stateAbbSort <- sort(aegypti$state_code[!duplicated(aegypti$state_code)]) #sorted vector of abbreviation

zikaModel <- data.frame()

plot.new()
plot.window(
  xlim = c(0, max(AC_rips$death)),
  ylim = c(0, max(AC_rips$death)),
  asp = 1
)
axis(1L)
axis(2L)
abline(a = 0, b = 1)
points(AC_rips[, c("birth", "death")], pch = 16L)

topologicalFeatures <- function(state_rips){
  H0_H1 <- ifelse(state_rips$dimension == 0, 0, 1)
  numH0 <- length(which(H0_H1 == 0))
  numH1 <- length(which(H0_H1 == 1))
  if(numH1 != 0){
    maxPerst <- max(state_rips$persistence[(numH0 + 1) : (numH0 + numH1)])
    medianPerst <- median(state_rips$persistence[(numH0 + 1) : (numH0 + numH1)])
    
  }else{
    maxPerst <- 0
    medianPerst <- 0
  }
  return(c(numH0, numH1, medianPerst, maxPerst))
}

for(val in stateAbbSort){
  state_coord <- aegypti[aegypti$state_code == val, c("x", "y"), drop = FALSE]
  state_rips <- vietoris_rips(state_coord) ##filtration
  state_rips$persistence<-(state_rips$death - state_rips$birth)
  
  features <- topologicalFeatures(state_rips)
  zikaModel <- rbind(zikaModel, c(features[1], features[2], features[3], features[4]))
}
colnames(zikaModel) <-  c("H0", "H1", "H1_Median","H1ML")
rownames(zikaModel) <- stateAbbSort

zikaModel <- merge(zikaModel, zikaPredict, by = "row.names")[,-1]
rownames(zikaModel) <- stateAbbSort

#-----Plot of Cases Across States-----#
plot(zikaModel$CaseNum, xaxt = "n", xlab = "States", ylab = "Number of Cases")
axis(1, at = 1:27, labels = stateAbbSort)

hist(log(zikaModel$CaseNum), main = "Distribution of Cases with Log Transformation")

#-----Linear Regression zikaModel----#
#We use the following plot to check if there is any correlation between predictors. 
pairs(zikaModel)

#-----------Model 1: without topological features-------------#
#First, we fit population, temperature, precipitation, and H0 to linear model with log transformation. 
fit.1 <- lm(log(CaseNum) ~ POP + TEMP + Precip + H0, data = zikaModel);summary(fit.1)
plot(fit.1, which = 1)
plot(fit.1, which = 2)
#Based on the summary output, F = 4.059 and p  = 0.01363. Coefficients of these predictors are not all 0.
#Checking assumptions of error, we plot residuals vs. fitted values and theoretical quantiles vs. standardized residuals. We observe that residuals for each fitted values are evenly distributed above and below the mean 0, and the relationship between residuals and fitted values seem to have a linear relation. 
#In the qauntile-quantile plo,t standardized residuals closely line up with normal values, with a few outliers. Therefore, the normal distribution assumption on error term is satisfied. 

#We use the liklihood ratio test to assess if we should stick with simpler model. 
fit.nested <- lm(log(CaseNum) ~ POP + TEMP + H0, data = zikaModel)
lmtest::lrtest(fit.nested, fit.1)
fit.1 <- lm(log(CaseNum) ~ POP + TEMP + H0, data = zikaModel);summary(fit.1)
#In the liklihood ratio test, p = 0.8135, which fails to reject the nested model. This result suggests that the simpler model, without "Precip", is better. 

#------------Model 2 : with topological features---------------#
#Then, we want to create a linear model with topological features as predictors to compare with "fit.1". Based on the scatterplot matrices above, we observe that H0 and H1 are highly correlated. We are concerned about effects of the collinearity on the linear regression, so we wonder the proportion of observed variation in "H1" explained by other predictors.
colinearModel <- lm(log(CaseNum) ~ . -H1, data = zikaModel);summary(colinearModel)$r.squared
#We see that other predictors explain 50.83% of variations in H1, thus we decided to keep H1 in our model before proceeding to variable selection. 

#Let's add H1, H0, H1ML, and H1_Median to our regression model. Similarly, error assumptions are satisfied by using two plots below.
fit.2 <- lm(log(CaseNum) ~ POP + TEMP + Precip + H0 + H1 + H1ML + H1_Median, data = zikaModel);summary(fit.2)
plot(fit.2, which = 1)
plot(fit.2, which = 2)

#Based on the summary output, F = 3.268 and p = 0.01878. We know that coefficients of these predictors are not all 0.
#However, We might overfit the data as many predictors are not statistically significant. Thus, we want to leave out variables that are unnecessary. 

#Using likelihood ratio test, we fail to reject the nested model when H1_Median is absent (p = 0.5024). Thus, we can leave out H1_Median in our model.  
fit.test0 <- lm(log(CaseNum) ~ POP + TEMP + Precip + H0 + H1 + H1ML, data = zikaModel)
lmtest::lrtest(fit.test0, fit.2)

#Also, we want to investigate if it is necessary to include TEMP and Precip in our model, as their coefficients are not statistically significant. 
fit.test1 <- lm(log(CaseNum) ~ POP + H0 + H1 + H1ML, data = zikaModel);summary(fit.test1)
lmtest::lrtest(fit.test1, fit.test0)
#In the likelihood ration test , we fail to reject the null (p = 0.2702). We choose the nested model that leaves out TEMP and Precip.

#We observe that coefficients of both H0 and H1 are not statistically significant. Also, it is concerning that the coefficients of H0 is negative, which means more mosquitoes are likely to decrease the number of cases. 
#We decide to test the model without H1, as it seems to correlated with H0.
fit.test2 <- lm(log(CaseNum) ~ POP + H0 + H1ML, data = zikaModel);summary(fit.test2)
lmtest::lrtest(fit.test2, fit.test1)
#The likelihood ratio test suggests that we need to choose the simpler model without H1 (p = 0.2364).

fit.2 <- lm(log(CaseNum) ~ POP + H0 + H1ML, data = zikaModel)
summary(fit.2)
summary(fit.1)

plot(x = exp(predict(fit.2)), y = zikaModel$CaseNum, xlab = "Number of Cases (Predicted)", ylab = "Number of Cases (Actual)")
abline(a = 0, b = 1)

plot(x = exp(predict(fit.1)), y = zikaModel$CaseNum, xlab = "Number of Cases (Predicted)", ylab = "Number of Cases (Actual)")
abline(a = 0, b = 1)

Example_States <- c("PE", "CE", "GO")
Coord_Example <- data.frame()
for(val in Example_States){
  targetRow <- zikaModel[rownames(zikaModel) == val,]
  Coord_Example <- rbind(Coord_Example, c(targetRow$H0, targetRow$H1ML, targetRow$CaseNum))
}
rownames(Coord_Example) <- Example_States
colnames(Coord_Example) <- c("H0", "H1ML", "CaseNum");Coord_Example


PE_coord <- aegypti[aegypti$state_code == "PE", c("x", "y"), drop = FALSE]
PE_rips <- vietoris_rips(PE_coord)
plot(x = PE_coord$x, y = PE_coord$y, col = "green", xlim = c(-42, -32))

CE_coord <- aegypti[aegypti$state_code == "CE", c("x", "y"), drop = FALSE]
CE_rips <- vietoris_rips(CE_coord)
plot(x = CE_coord$x, y = CE_coord$y, col = "blue", xlim = c(-43, -33))

GO_coord <- aegypti[aegypti$state_code == "GO", c("x", "y"), drop = FALSE]
GO_rips <- vietoris_rips(GO_coord) 
plot(x = GO_coord$x, y = GO_coord$y, col = "red", xlim = c(-54, - 44))

