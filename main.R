library(xtable)
library(generalhoslem)
library(ggplot2)
library(MCMCpack)
library(MASS)
library(dplyr)

setwd("C:/Users/duviv/Documents/University/KUL/S2/Generalized Linear Models/Individual Homework")
df <- read.csv("student-mat.csv")
df1 <- df[,c("Score", "sex", "age", "Pstatus", "Medu", "famsup", "paid", 
            "higher", "internet", "romantic")]
str(df1) #medu needs to be transformed as a factor
df1$Medu <- as.factor(df1$Medu)
str(df1$Medu)
df1$Score <- as.factor(ifelse(df1$Score == "Pass", df1$Score <- 1, df1$Score <- 0))
df1$higher <- as.factor(ifelse(df1$higher == "yes", df1$higher <- 1, df1$higher <- 0))

sum(is.na(df[,1:10])) #Not one single NaN

# 0 - Descriptive statistics ----

xtabs(~ Score + sex, data = df)
xtabs(~ Score + age, data = df)
xtabs(~ Score + Pstatus, data = df)
xtabs(~ Score + Medu, data = df)
xtabs(~ Score + famsup, data = df)
xtabs(~ Score + paid, data = df)
xtabs(~ Score + higher, data = df)
xtabs(~ Score + internet, data = df)
xtabs(~ Score + romantic, data = df)

xtable(df1)


####
# ANALYSIS 1 ----
####

# 1 - Explain Score with age, Medu, paid, higher, and romantic ----
# 1a - Frequentist manner ----

#fit the full model
df.freq.fit1 <- glm(data = df1, Score ~ age + Medu + paid + higher + romantic , family = binomial("logit"))
summary(df.freq.fit1) # Age is a continuous variable so Bernouilli assumption is needed here
str(df.freq.fit1)
# (1) which covariates to include in the model ----

#variable selection, based on the Rao test
step(df.freq.fit1, test = "Rao") #Suggests to keep the intercept, age, and higher
df.freq.fit2 <- glm(data = df1, Score ~ age + higher, family = binomial("logit"))
summary(df.freq.fit2)

coef <- coefficients(df.freq.fit2)
n1 <- c(1, 16, 1)
n2 <- c(1, 20, 0)

1/(1 + exp(-sum(n1*coef)))
1/(1 + exp(-sum(n2*coef)))

# (2) model adequacy, ----

# 2.1 Goodness of Fit ----
# We have the continuous variable "age" so should use the Homer-Lemeshow test
# install.packages("generalhoslem")

logitgof(obs = df1$Score, exp = fitted(df.freq.fit2), g=10, ord = FALSE)
# hosmerlem = function(y, yhat, g=10) # Keep g=10 as this is the most commonly used number
#   {
#   fcutyhat = cut(yhat, breaks = quantile(yhat, probs=seq(0,1, 1/g)), include.lowest = TRUE)
#   obs = xtabs(cbind(1 - y, y) ~ cutyhat)
#   expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
#   chisq = sum((obs - expect)^2/expect)
#   P = 1 - pchisq(chisq, g - 2)
# 
#   return(list(chisq=chisq,p.value=P))
# }
# hosmerlem(y=df1$Score, yhat=fitted(df.freq.fit3))

# 2.2 Quality of Prediction----

## Now calculate the overall "Pseudo R-squared" and its p-value
ll.null <- df.freq.fit2$null.deviance/-2
ll.proposed <- df.freq.fit2$deviance/-2

## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
pseudo.R2 <- (ll.null - ll.proposed) / ll.null

## The p-value for the R^2, test if this model is better than a model with only an intercept
p.val <- 1 - pchisq(2*(ll.proposed - ll.null), df=(length(df.freq.fit2$coefficients)-1))
# Reject the fact that both models are similar

## Display results
print(matrix(data = c("Pseudo R²", "Overall p-value", pseudo.R2, p.val), nrow = 2, ncol = 2, byrow = TRUE))

# Calculate Tjur's proposal
tjur <- cbind(fitted(df.freq.fit2), df1$Score)
tjur0 <- subset(tjur[,1], tjur[,2]==1)
tjur1 <- subset(tjur[,1], tjur[,2]==2)
mean0 <- mean(tjur0)
mean1 <- mean(tjur1)
D <- mean1 - mean0

# 2.3 Residual diagnostic----
# Outlier diagnostic
r.dev <- residuals(df.freq.fit2, type = "deviance")
plot(df1$age, r.dev, xlab="Age (years)", ylab="Deviance residual", cex.lab=1.5, cex.axis=1.3, 
     main = "Residual plot of deviance residual versus age")
loess.dev <- loess(r.dev ~ df1$age)
lo.pred <- predict(loess.dev, se=T)
orderage <- order(df1$age)
lines(df1$age[orderage],lo.pred$fit[orderage], col="blue", lwd=3)
lines(df1$age[orderage],lo.pred$fit[orderage] + 2*lo.pred$s[orderage], lty=2, col="red")
lines(df1$age[orderage],lo.pred$fit[orderage] - 2*lo.pred$s[orderage], lty=2, col="red")

# Influential observations diagnostic
N <- length(df1$Score)
id <- 1:N
# lrmod2 <- glm(bpd~birthwt,family=binomial,data=bw2)

# Global influence plots
hat.df1 <- hatvalues(df.freq.fit2)
rstudent.df1 <- rstudent(df.freq.fit2)
plot(hat.df1, rstudent.df1, main = "Scatter plot of hat-values versus studentized residuals")
dffits.df1 <- dffits(df.freq.fit2)
plot(id, dffits.df1, type="l", main = "Index plot of Dffits")
cov.df1 <- covratio(df.freq.fit2)
plot(id,cov.df1, type="l", main = "Index plot of Covariance ratio")
cook.df1 <- cooks.distance(df.freq.fit2)
plot(id, cook.df1, type="l", main = "Index plot of Cook's distance")

influencers <- subset(cook.df1, cook.df1> 0.01)
influencers
influencers.id <- as.integer(names(influencers)) # Do not take them out because we do not have justification to do so
df1[influencers.id, 3]

# (3) and predict the response using the final model ----

## Lastly, let's  see what this logistic regression predicts, given
## that a patient is either female or male (and no other data about them).
predicted.data <- data.frame(
  probability.to.pass = df.freq.fit2$fitted.values,
  Score = df1$Score, higher = df1$higher, age = df1$age)

ggplot(data = predicted.data, aes(x=Score, y = probability.to.pass)) +
  geom_point(aes(color = Score), size = 2) +
  xlab("Score") +
  ylab("Predicted probability to pass")

## now we can plot the data
predicted.data <- data.frame(
  probability.to.pass = df.freq.fit2$fitted.values,
  Score = df1$Score, higher = df1$higher, age = df1$age) 

predicted.data <- predicted.data[order(predicted.data$probability.to.pass, decreasing=FALSE),]
predicted.data$rank <- 1:nrow(predicted.data)

## Lastly, we can plot the predicted probabilities for each sample having
## heart disease and color by whether or not they actually had heart disease
ggplot(data = predicted.data, aes(x=rank, y=probability.to.pass)) +
  geom_point(aes(color=Score), alpha=1, stroke=1) +
  xlab("Index") +
  ylab("Predicted probability to pass")

ggsave("passing_probabilities3.jpg")

# 1b - Bayesian manner ----
# You are only asked to fit the finally chosen frequentist model with a Bayesian program
score.bayes <- MCMClogit(Score ~ age + higher,family=binomial, data=df1)
summary(score.bayes)

# 1c - Interpret the results ----
# De la rédaction de ouf

####
# ANALYSIS 2 ----
####
setwd("C:/Users/duviv/Documents/University/KUL/S2/Generalized Linear Models/Individual Homework")

df2 <- read.csv("student-mat.csv")
df2 <- df2[,c("absences", "sex", "Pstatus", "famsup", "internet")]
str(df2)
df2$sex <- ifelse(df2$sex=="M", 1, 0)
df2$Pstatus <- ifelse(df2$Pstatus=="T", 1, 0)
df2$famsup <- ifelse(df2$famsup=="yes", 1, 0)
df2$internet <- ifelse(df2$internet=="yes", 1, 0)


# 2a - Relate the number of school absences (numeric) to sex, Pstatus, famsup and internet----
fit.poiss <- glm(data = df2, absences ~ ., family = poisson(link="log"))
summary(fit.poiss)
xtable(summary(fit.poiss))

#we compare with NB regression to see if there is any improvement
fit.nb <- glm.nb(data = df2, absences ~ .)
summary(fit.nb) #AIC is 2188 instead of 4092 so we keep the nb model
xtable(summary(fit.nb))

step(fit.nb, test = "Rao") #Keep intercept Pstatus and internet
xtable(step(fit.nb, test = "Rao"))

# 2c - Check in a frequentist manner: ----
  # (1) which covariates to include in the model, ----
step(fit.nb, test = "Rao")
fit.nb.2 <- glm.nb(data = df2, absences ~ Pstatus + internet)
summary(fit.nb.2)

#test students
1/(1 + exp(-(1.7984)))


  # (2) model adequacy ----
# 2.1 - Goodness of fit ----

dev <- summary(fit.nb.2)$deviance # Almost no difference of deviance between the full model and the selected model
df <- summary(fit.nb.2)$df.residual
1-pchisq(dev,df) #which means our model is better than the null model

# 2.2 - Adequacy of prediction ----

## Now calculate the overall "Pseudo R-squared" and its p-value
ll.null <- fit.nb.2$null.deviance/-2
ll.proposed <- fit.nb.2$deviance/-2

## McFadden's Pseudo R^2
pseudo.R2 <- (ll.null - ll.proposed) / ll.null

## The p-value for the R^2, test if this model is better than a model with only an intercept
p.val <- 1 - pchisq(2*(ll.proposed - ll.null), df=(length(fit.nb.2$coefficients)-1))
# Reject the fact that both models are similar

## Display results
print(matrix(data = c("Pseudo R²", "Overall p-value", pseudo.R2, p.val), nrow = 2, ncol = 2, byrow = TRUE))

# 2.3 - Residual diagnostic ----

# Outlier diagnostic
r.dev <- residuals(fit.nb.2, type = "deviance")
plot(df2$absences, r.dev, xlab="Days of absence", ylab="Deviance residual", cex.lab=1.5, cex.axis=1.3, 
     main = "Residual plot of deviance residual versus age")
loess.dev <- loess(r.dev ~ df2$absences)
lo.pred <- predict(loess.dev, se=T)
orderage <- order(df2$absences)
lines(df2$absences[orderage],lo.pred$fit[orderage], col="blue", lwd=3)
lines(df2$absences[orderage],lo.pred$fit[orderage] + 2*lo.pred$s[orderage], lty=2, col="red")
lines(df2$absences[orderage],lo.pred$fit[orderage] - 2*lo.pred$s[orderage], lty=2, col="red")

## Influential observations diagnostic
N <- length(df2$absences)
id <- 1:N

# Global influence plots
hat.df2 <- hatvalues(fit.nb.2)
rstudent.df2 <- rstudent(fit.nb.2)
plot(hat.df2, rstudent.df2, main = "Scatter plot of hat-values versus studentized residuals")
dffits.df2 <- dffits(fit.nb.2)
plot(id, dffits.df2, type="l", main = "Index plot of Dffits")
cov.df2 <- covratio(fit.nb.2)
plot(id,cov.df2, type="l", main = "Index plot of Covariance ratio")
cook.df2 <- cooks.distance(fit.nb.2)
plot(id, cook.df2, type="l", main = "Index plot of Cook's distance")

# influencers <- subset(cook.df2, cook.df2> 0.01)
# influencers
# id.infl <- c(as.integer(names(influencers)))
# df2[id.infl,1]
# 
# influencers <- data.frame(id = as.integer(names(influencers)), 
#                       cook = as.integer(subset(cook.df2, cook.df2> 0.01)),
#                       days = df2[id.infl,1])
#                       # nrow = length(influencers)
#                       # ncol = 3, 
#                       # byrow = FALSE)
# 
# influencers.id <- as.integer(names(influencers)) # Do not take them out because we do not have justification to do so
# df2$absences[influencers.id]

  # (3) predict the response using the final model ----

xtable(table(df2$absences))
counts <- table(df2$absences)
barplot(counts,
        xlab="Days of absence",
        ylim = c(0, 120))

predicted.data <- data.frame(
  days.predicted = fit.nb.2$fitted.values,
  Absences = df2$absences, 
  Parental_Status = df2$Pstatus, 
  Internet = df2$internet)

ggplot(data = predicted.data, aes(x=Absences, y = days.predicted)) +
  geom_point(aes(color = Absences), size = 2) +
  xlab("Observed days of absence") +
  ylab("Predicted days of absence") +
  geom_abline(color = "red", slope = 1, intercept = 0)

# 2d - For the Bayesian analysis you are only asked to fitt the finally chosen frequentist model with a Bayesian program----
fit.nb.bayes <- MCMCpoisson(absences ~ Pstatus + internet, data = df2)
summary(fit.nb.bayes)