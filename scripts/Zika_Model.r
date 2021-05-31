devtools::load_all()
data(aegypti)
data(zikaPredict)

stateAbbSort <- sort(aegypti$state_code[!duplicated(aegypti$state_code)]) #sorted vector of abbreviation
stateAbbSort <- stateAbbSort[stateAbbSort != "DF"]

zikaModel <- data.frame()

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
axis(1, at = 1:26, labels = stateAbbSort)

hist(log(zikaModel$CaseNum), main = "Distribution of Cases with Log Transformation")

#-----Linear Regression zikaModel----#
pairs(zikaModel)

#-----------Model 1: without topological features-------------#
fit.1 <- lm(log(CaseNum) ~ POP + TEMP + H0 + Precip, data = zikaModel);summary(fit.1)
plot(fit.1, which = 1)
plot(fit.1, which = 2)
anova(fit.1)
#Based on the summary output, F = 4.059 and p  = 0.01363. Coefficients of these predictors are not all 0.
#Based on the ANOVA table, we can leave out Precip as a predictor (p = 0.6248)
fit.1 <- lm(log(CaseNum) ~ POP + TEMP + H0, data = zikaModel);summary(fit.1)

#------------Model 2 : with topological features---------------#
fit.2 <- lm(log(CaseNum) ~ POP + TEMP + H0 + H1ML + H1 + H1_Median, data = zikaModel);summary(fit.2)
plot(fit.2, which = 1)
plot(fit.2, which = 2)
anova(fit.2)
#Based on the summary output, F = 6.129 and p = 0.002366. Coefficients of these predictors are not all 0.
#We might overfit the data
#Based on the ANOVA table, we can leave out H1_Median as a predictor (p = 0.3262)

fit.2 <- lm(log(CaseNum) ~ POP + TEMP + H0 + H1ML + H1, data = zikaModel);summary(fit.2)
fit.test0 <- lm(log(CaseNum) ~ POP + H0 + H1ML + H1, data = zikaModel);summary(fit.test0)
lmtest::lrtest(fit.test0, fit.2)
#Choose nested model (without TEMP)
fit.test1 <- lm(log(CaseNum) ~ POP + H0 + H1ML, data = zikaModel);summary(fit.test1)
lmtest::lrtest(fit.test1, fit.test0)
#Choose nested model (without H1)

fit.2 <- lm(log(CaseNum) ~ POP + H0 + H1ML, data = zikaModel);summary(fit.2)


plot(x = exp(predict(fit.2)), y = zikaModel$CaseNum, xlab = "Number of Cases (Predicted)", ylab = "Number of Cases (Actual)")
abline(a = 0, b = 1)

plot(x = exp(predict(fit.1)), y = zikaModel$CaseNum, xlab = "Number of Cases (Predicted)", ylab = "Number of Cases (Actual)")
abline(a = 0, b = 1)


#----------------------Example: Topological Features in State "GO"----------------------#
state_coord <- aegypti[aegypti$state_code == "GO", c("x", "y"), drop = FALSE]
state_rips <- vietoris_rips(state_coord) ##filtration
plot.new()
plot.window(
  xlim = c(0, max(state_rips$death)),
  ylim = c(0, max(state_rips$death)),
  asp = 1
)
axis(1L)
axis(2L)
abline(a = 0, b = 1)
points(state_rips[, c("birth", "death")], pch = 16L)

#Plot of Aedes Aegypti Distribution Within a State
plot(x = state_coord$x, y = state_coord$y)
