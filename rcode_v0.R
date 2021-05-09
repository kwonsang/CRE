#######################################
#### Discovery step of the CRE method
#######################################
## Read the dataset
dataset = read.csv("simulated_dataset.csv")

y = dataset[,1] # outcome
z = dataset[,2] # treatment
X = dataset[,3:12] # 10 covariates

## There are 10 covariates (x1, x2, ..., x10).
## The correlation between x1, x2, ..., x10 is 0.1.

## The true propensity score model log(p/(1-p)) = -1 + x1 - x2 + x3.
## treatment z is generated based on a Bernoulli distribution with p.

## The outcome follows a normal distribution. 
## There is effect modification by x1, x2, x3
## The potential outcome y(0) ~ N(0,1)
## y(1) = y(0) + tau
## If x1=0 & x2=0, tau=2
## If x1=1 & x3=1, tau=-2
## Otherwise, tau=0

#######################################
source("causal_rulefit_functions.R")

## 1. Rule generation
### 1.1. obtain an estimate for ITE
#### Use BCF
library(bcf)

X = as.matrix(X)
propscore.model = glm(z ~ X, family = binomial)
logit.ps = predict(propscore.model)
est.ps = exp(logit.ps)/(1+exp(logit.ps))

bcf.model = bcf(y, z, X, X, est.ps, nburn = 100, nsim = 1000)
tau.est = colMeans(bcf.model$tau)

### 1.2. fit tree ensemble
### 1.3. extract decision rules
library(randomForest)
library(gbm)
library(inTrees)
library(xgboost)

S=20
L=5
ntree=200
### rule extraction & regularization!
y.temp = tau.est
muy = mean(y.temp)
sdy = sd(y.temp)
yz = (y.temp-muy)/sdy # standardize the treatment effect

rulesf = genrulesRF(X, yz, nt=ntree, S=S, L=L) # rule set from RF
rules.gbm = genrulesGBM(X, yz, nt=ntree, S=S, L=L) # rule set from GBM

# combine the extracted rule sets
rulesf = c( rulesf, rules.gbm)
minsup=0.025 # minimum support
dt = createX.take1(X= X, rules = rulesf, t = minsup)
#dt = createX.disjoint(X= X, rules = rulesf, t = minsup)
Xr = dt[[1]]
rulesFin = dt[[2]]

## 2. Rule-regularization 
### 2.1. Generate tilde{X}
std.Xr = Xr
mur = apply(Xr, 2, mean)
sdr = apply(Xr, 2, sd)
for(l in 1:dim(Xr)[2]){
  std.Xr[,l] = (Xr[,l]-mur[l])/sdr[l]
}

### 2.2. Penalized regression
#### 2.2.1. LASSO


library(glmnet)
lambda = 10^seq(10, -2, length = 100)
lasso.mod = glmnet(std.Xr, yz, alpha = 1, lambda = lambda, intercept = FALSE)
cv.lasso = cv.glmnet(std.Xr, yz, alpha = 1, intercept = FALSE)
bestlam = cv.lasso$lambda.min

aa=coef(cv.lasso, s=cv.lasso$lambda.1se)
index.aa = which(aa[-1,1]!=0)

rule.LASSO = data.frame(rules = rulesFin[index.aa], val = aa[index.aa+1, 1])
rule.LASSO[order(-rule.LASSO[,2]), ]


#### 2.2.2. Stability selection 

library(stabs)
stab.mod = stabsel(std.Xr, yz, fitfun = glmnet.lasso, cutoff = 0.8, PFER=1, args.fitfun = list(type = "conservative"))
plot(stab.mod, main = "Lasso")
rule.stab = rulesFin[stab.mod$selected]
rule.stab 


