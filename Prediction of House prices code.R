install.packages("leaps")
install.packages("glmnet")
install.packages("corrplot")
install.packages("lasso2")
library(lasso2)
library(glmnet)
library(MASS)
library(GGally)
library(scales)
library(dplyr)
library(leaps)
library(corrplot)

?Boston
Boston1<-Boston
View(Boston1)
str(Boston1)
summary(Boston1)
unique(Boston1$rad)
#Correlation
ggpairs(Boston1[-c(4)])

as.matrix()
corrplot(cor(Boston1[-c(4)]), type="upper", order="hclust", col=c("black", "white"),bg="lightblue")
corrplot(cor(Boston1[-c(4)]), method="circle")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
par(mfrow=c(1,1))
# matrix of the p-value of the correlation
p.mat <- cor.mtest(Boston1[-c(4)])
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(Boston1[-c(4)]), method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

#as Factors
Boston1$chas<-as.factor(chas)
str(Boston1)

par(mfrow=c(1,1))
#Boxplots for outliers
boxplot(Boston,col = c("skyblue"), main ="Box Plot for all variables")

#Standardize
Boston<-scale(Boston1[-c(4,14)],center = FALSE, scale = apply(Boston1[-c(4,14)],2, sd, na.rm = TRUE))

Boston<-data.frame(Boston)
Boston$chas<-Boston1$chas
Boston$medv<-Boston1$medv
View(Boston)
attach(Boston)


# 70 % training data
sample_index <- sample(nrow(Boston),nrow(Boston)*0.75)
Boston_train <- Boston[sample_index,]
Boston_test <- Boston[-sample_index,]


#Initial Models - 13 variables (all)

model_1 <- lm(medv~., data=Boston_train)
model_1_summ<-summary(model_1)
model_1$call
fit.mse <- (model_1_summ$sigma)^2
fit.rsq <- model_1_summ$r.squared
fit.arsq <- model_1_summ$adj.r.squared
test.pred.fit <- predict(model_1, newdata=Boston_test) 
fit.mpse <- mean((Boston_test$medv-test.pred.fit)^2)
fit.aic <- AIC(model_1)
fit.bic <- BIC(model_1)

stats.fit <- c("FULL", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)

comparison_table <- c("model type", "MSE", "R-Squared", "Adjusted R-Squared", "Test MSPE", "AIC", "BIC")
data.frame(cbind(comparison_table, stats.fit))

#P-values for the variables indus and age is insignificant so eliminated

model_2 <- lm(medv~crim+zn+chas+nox+rm+dis+rad+tax+ptratio+black+lstat, data=Boston_train)
model_2_summ<-summary(model_2)


fit.mse <- (model_2_summ$sigma)^2
fit.rsq <- model_2_summ$r.squared
fit.arsq <- model_2_summ$adj.r.squared
test.pred.fit <- predict(model_2, newdata=Boston_test) 
fit.mpse <- mean((Boston_test$medv-test.pred.fit)^2)
fit.aic <- AIC(model_2)
fit.bic <- BIC(model_2)

stats.p.fit <- c("P-Value Model", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)

data.frame(cbind(comparison_table, stats.fit,stats.p.fit))



##Feature selection
#Best subset model

#regsubsets only takes data frame as input
subset_result <- regsubsets(medv~.,data=Boston_train, nbest=2, nvmax = 14)
summary(subset_result)
plot(subset_result, scale="bic")

model_3 <- lm(medv~crim+zn+chas+nox+rm+dis+rad+ptratio+black+lstat, data=Boston_train)
model_3_summ<-summary(model_3)


fit.mse <- (model_3_summ$sigma)^2
fit.rsq <- model_3_summ$r.squared
fit.arsq <- model_3_summ$adj.r.squared
test.pred.fit <- predict(model_3, newdata=Boston_test) 
fit.mpse <- mean((Boston_test$medv-test.pred.fit)^2)
fit.aic <- AIC(model_3)
fit.bic <- BIC(model_3)

stats.ss.fit <- c("Best subsets Model", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)

data.frame(cbind(comparison_table, stats.fit,stats.p.fit,stats.ss.fit))


#If there are n independent variables,
#the number of possible nonempty subsets is 2^n - 1. 
#If you try a best subset regression with more than 50 variables, 
#you might need to wait for your entire life to get the result.

#Stepwise selection
nullmodel=lm(medv~1, data=Boston_train)
fullmodel=lm(medv~., data=Boston_train)

#Forward
model_step_f <- step(nullmodel, scope=list(lower=nullmodel, upper=fullmodel), direction='forward')
model_3_summ<-summary(model_step_f)

stats.f.fit <- c("Forward selection", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)

data.frame(cbind(comparison_table, stats.fit,stats.p.fit,stats.ss.fit,stats.f.fit))


#Backward - same as p value and subset and both step
model_step_b <- step(fullmodel,direction='backward')
model_4_summ<-summary(model_step_b)

stats.b.fit <- c("Backward Elimination", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)
data.frame(cbind(comparison_table, stats.fit,stats.p.fit,stats.ss.fit,stats.f.fit,stats.b.fit))

#Both
model_step_s <- step(nullmodel, scope=list(lower=nullmodel, upper=fullmodel), direction='both')
#same as backward

stats.step.fit <- c("Stepwise", fit.mse, fit.rsq, fit.arsq, fit.mpse, fit.aic, fit.bic)
data.frame(cbind(comparison_table, stats.fit,stats.p.fit,stats.ss.fit,stats.f.fit,stats.b.fit,stats.step.fit))


## Lasso
lasso_fit = glmnet(x = as.matrix(Boston_train[, -c(which(colnames(Boston_train)=='medv'))]), y = Boston_train$medv, alpha = 1)
plot(lasso_fit, xvar = "lambda")
#lambda = 0.5
coef(lasso_fit,s=0.5)

#lambda = 1
coef(lasso_fit,s=1)

#use 5-fold cross validation to pick lambda
cv_lasso_fit = cv.glmnet(x = data.matrix(Boston_train[, -c(which(colnames(Boston_train)=='medv'))]), y = Boston_train$medv, alpha = 1, nfolds = 5)
plot(cv_lasso_fit)

cv_lasso_fit$lambda.min
cv_lasso_fit$lambda.1se

#lambda = 0.36
coef(lasso_fit,s=cv_lasso_fit$lambda.min)

mse = function(x,y) { mean((x-y)^2)}

#MSE
lasso.train = predict(cv_lasso_fit, s = cv_lasso_fit$lambda.min, newx = data.matrix(Boston_train[, -c(which(colnames(Boston_train)=='medv'))])) 
mse.1=mse(lasso.train,Boston_train$medv)

#Prediction error
lasso.pred = predict(cv_lasso_fit, s = cv_lasso_fit$lambda.min, newx = data.matrix(Boston_test[, -c(which(colnames(Boston_test)=='medv'))])) 
mspe.1 = mse(lasso.pred,Boston_test$medv)

#R sqr

rsq <- lasso_fit$dev.ratio[which(cv_lasso_fit$glmnet.fit$lambda==cv_lasso_fit$lambda.min)]
n=506
arsq <- 1-((1-rsq)*(n-1)/(n-13-1))
#

stats.l.fit <- c("Lasso Model", mse.1, rsq, arsq, mspe.1, "N/A", "N/A")

data.frame(cbind(comparison_table, stats.fit,stats.p.fit,stats.ss.fit,stats.f.fit,stats.b.fit,stats.step.fit,stats.l.fit))

par(mfrow=c(2,2))
plot(model_step_s)




