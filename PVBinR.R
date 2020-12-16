# Data, 2003a Kosinski & Barnhart
# X1 = Gender: 1 = male, 0 = female.
# X2 = Stress mode: 1 = dipyridamole, 0 = exercise.
# X3 = Age: 1 = age ≥ 60 years, 0 = age < 60 years.
# T = SPECT thallium test: 1 = positive, 0 = negative.
# D = CAD (disease) status: 1 = yes, 0 = no; NA = not verified.
cad = read.csv("2003kosinski_cad.csv")
str(cad)

# Verification status
cad$R = 1  # R = Verified: 1 = yes, 0 = no
cad[is.na(cad$D), "R"] = 0

# Percent missing
mean(is.na(cad$D))*100

# Library
library(prediction)
library(mice)

# Session Info
sessionInfo()

# Complete-Case Analysis
tbl = table(cad[, c("T", "D")]); tbl
sn0 = tbl[2,2]/sum(tbl[,2]); sn0
sp0 = tbl[1,1]/sum(tbl[,1]); sp0

# Recode NA to -1 for ease of calculation for B&G
cad_ = cad
cad_[is.na(cad_$D), "D"] = -1  # recode NA as -1, this will be -1*0=0 for unverified. Else NA*0=NA.

# B&G, Count Method in Begg & Greenes 1983
tbl_ = table(cad_[, c("T", "D")]); addmargins(tbl_)  # -1 column -- unverified cases
upper = (sum(cad_$T)/nrow(cad_)) * (sum(cad_$D*cad_$T*cad_$R)/sum(cad_$T*cad_$R))
lower = (sum(1-cad_$T)/nrow(cad_)) * (sum(cad_$D*(1-cad_$T)*cad_$R)/sum((1-cad_$T)*cad_$R))
sn1 = upper / (upper+lower); sn1
upper1 = (sum(1-cad_$T)/nrow(cad_)) * (sum((1-cad_$D)*(1-cad_$T)*cad_$R)/sum((1-cad_$T)*cad_$R))
lower1 = (sum(cad_$T)/nrow(cad_)) * (sum((1-cad_$D)*cad_$T*cad_$R)/sum(cad_$T*cad_$R))
sp1 = upper1 / (upper1+lower1); sp1

# B&G, Regression Method in Alonzo 2005
model1 = glm(D ~ T, data = cad, family = "binomial")  # Modelled using complete data
summary(model1)
preds = prediction(model1, cad_)$fitted  # Predicted on incomplete data
sn2 = sum(cad_$T*preds) / sum(preds); sn2            # P(T=1|D=1)
sp2 = sum((1-cad_$T)*(1-preds)) / sum(1-preds); sp2  # P(T=0|D=0)

# B&G, Regression Method in Alonzo 2005, with one covariate
model1x = glm(D ~ T*X1, data = cad, family = "binomial")  # Modelled using complete data
summary(model1x)
preds = prediction(model1x, cad_)$fitted  # Predicted on incomplete data
sn2x = sum(cad_$T*preds) / sum(preds); sn2x            # P(T=1|D=1,X1)
sp2x = sum((1-cad_$T)*(1-preds)) / sum(1-preds); sp2x  # P(T=0|D=0,X1)

# Multiple Imputation, Logistic Regression
# Modified codes from deGroot et al 2008
# Logistic regression
m = 85  # number of imputed data sets
seednum = 12345
cad1 = cad[, c("T","D")]  # select only T & D, else MI use all variables for imputation
cad1$D = as.factor(cad1$D); data_impute = mice(cad1, m = m, method = "logreg", seed = seednum)  # logistic regression
# results dependent on seed & method, seems pmm less variant than logreg
data_impute$method
cad_mi = complete(data_impute, action = "repeated")  # MI data
str(cad_mi)
# loop to calculate sn & spe
sn3m = sp3m = rep(NA, m)
for (i in 1:m) {
  tbl = table(T=cad_mi[,1],D=cad_mi[,m+1])
  sn3m[i] = tbl[2,2]/sum(tbl[,2])
  sp3m = tbl[1,1]/sum(tbl[,1])
}
sn3 = mean(sn3m); sn3
sp3 = mean(sp3m); sp3
# MI by logistic regression, one covariate
m = 85  # number of imputed data sets = % missing observations
seednum = 12345
cad1 = cad[, c("T","D", "X1")]  # select only T, D and X1, else MI use all variables for imputation
cad1$D = as.factor(cad1$D); data_impute = mice(cad1, m = m, method = "logreg", seed = seednum)  # logistic regression
data_impute$method
cad_mi = complete(data_impute, action = "repeated")  # MI data
# names(cad_mi)
# loop to calculate sn & spe
sn3m = sp3m = rep(NA, m)
for (i in 1:m) {
  tbl = table(T=cad_mi[,1],D=cad_mi[,m+1])
  sn3m[i] = tbl[2,2]/sum(tbl[,2])
  sp3m = tbl[1,1]/sum(tbl[,1])
}
sn3x = mean(sn3m); sn3x
sp3x = mean(sp3m); sp3x

# PMM
m = 85  # number of imputed data sets
seednum = 12345
cad1 = cad[, c("T","D")]  # select only T & D, else MI use all variables for imputation
cad1$D = as.factor(cad1$D); data_impute = mice(cad1, m = m, method = "pmm", seed = seednum)  # logistic regression
# results dependent on seed & method, seems pmm less variant than logreg
data_impute$method
cad_mi = complete(data_impute, action = "repeated")  # MI data
str(cad_mi)
# loop to calculate sn & spe
sn3m = sp3m = rep(NA, m)
for (i in 1:m) {
  tbl = table(T=cad_mi[,1],D=cad_mi[,m+1])
  sn3m[i] = tbl[2,2]/sum(tbl[,2])
  sp3m = tbl[1,1]/sum(tbl[,1])
}
sn3a = mean(sn3m); sn3a
sp3a = mean(sp3m); sp3a
# PMM, one covariate
m = 85  # number of imputed data sets = % missing observations
seednum = 12345
cad1 = cad[, c("T","D", "X1")]  # select only T, D and X1, else MI use all variables for imputation
cad1$D = as.factor(cad1$D); data_impute = mice(cad1, m = m, method = "pmm", seed = seednum)  # logistic regression
data_impute$method
cad_mi = complete(data_impute, action = "repeated")  # MI data
# names(cad_mi)
# loop to calculate sn & spe
sn3m = sp3m = rep(NA, m)
for (i in 1:m) {
  tbl = table(T=cad_mi[,1],D=cad_mi[,m+1])
  sn3m[i] = tbl[2,2]/sum(tbl[,2])
  sp3m = tbl[1,1]/sum(tbl[,1])
}
sn3ax = mean(sn3m); sn3ax
sp3ax = mean(sp3m); sp3ax

# EM Algorithm, Konsinski & Barnhart, 2003, MNAR
# pseudo-data
cad$R = 1  # R = Verified: 1 = yes, 0 = no
cad[is.na(cad$D), "R"] = 0
cad_pseudo = rbind(cad, cad[cad$R == 0, ])  # replicate U observation
str(cad_pseudo)
sum(is.na(cad_pseudo$D))  # 2217*2 = 4434
table(cad_pseudo$T)
# create 0, 1 for 2U rows
cad_pseudo[(sum(cad$R)+1):nrow(cad), "D"] = 0        # 1st U rows
cad_pseudo[(nrow(cad)+1):nrow(cad_pseudo), "D"] = 1  # 2nd U rows
# view
head(cad_pseudo)
head(cad_pseudo[(sum(cad$R)+1):nrow(cad),])
head(cad_pseudo[(nrow(cad)+1):nrow(cad_pseudo),])
n_u = nrow(cad_pseudo)  # total cases + U
# index k
index_1 = which(cad$R == 1)  # verified
index_2 = (sum(cad$R)+1):nrow(cad)  # unverified U
index_3 = (nrow(cad)+1):nrow(cad_pseudo)  # unverified 2U
# M-1 model, components: a = disease, b = diagnostic, c = missing data mechanism
# initialize values
weight_k = rep(1, nrow(cad_pseudo))  # init weight = 1 for all
coef_1 = vector("list", 0)
max_t = 500  # may increase
# EM Algorithm Iteration
for (t in 1:max_t) {
  # a -- P(D)
  model1a = glm(D ~ 1, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model1a = summary(model1a)
  coef_1a = su_model1a$coefficients
  fitted_pa = model1a$fitted.values
  # b -- P(T|D)
  model1b = glm(T ~ D, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model1b = summary(model1b)
  coef_1b = su_model1b$coefficients
  fitted_pb = model1b$fitted.values
  # c -- P(R|T)
  model1c = glm(R ~ T + D, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model1c = summary(model1c)
  coef_1c = su_model1c$coefficients
  fitted_pc = model1c$fitted.values
  # all
  coef_1 = append(coef_1, list(list(coef_1a, coef_1b, coef_1c)))
  # weight
  fitted_ps = list(fitted_pa, fitted_pb, fitted_pc)
  ys = list(model1a$y, model1b$y, model1c$y)
  p0 = (fitted_ps[[1]]^ys[[1]]) * ((1 - fitted_ps[[1]])^(1 - ys[[1]]))  # P(D)
  p1 = (fitted_ps[[2]]^ys[[2]]) * ((1 - fitted_ps[[2]])^(1 - ys[[2]]))  # P(T|D)
  p2 = (fitted_ps[[3]]^ys[[3]]) * ((1 - fitted_ps[[3]])^(1 - ys[[3]]))  # P(R|T)
  pk = p0 * p1 * p2  # P(R,T,D|X)
  weight_k[index_2] = pk[index_2] / (pk[index_2] + pk[index_3])
  weight_k[index_3] = 1 - weight_k[index_2]
  # next iteration
  t = t + 1
}
coef_1[[max_t]]
sn4 = predict(model1b, list(D=1), type="response"); sn4      # P(T=1|D=1)
sp4 = 1 - predict(model1b, list(D=0), type="response"); sp4  # P(T=0|D=0)

# EM Algorithm, Konsinski & Barnhart, 2003, MNAR
# M-2 model, components: a = disease, b = diagnostic, c = missing data mechanism
# initialize values
weight_k = rep(1, nrow(cad_pseudo))  # init weight = 1 for all
coef_2 = vector("list", 0)
max_t = 500  # may increase
# EM Algorithm Iteration
for (t in 1:max_t) {
  # a -- P(D)
  model2a = glm(D ~ X1, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model2a = summary(model2a)
  coef_2a = su_model2a$coefficients
  fitted_pa = model2a$fitted.values
  # b -- P(T|D)
  model2b = glm(T ~ D + X1, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model2b = summary(model2b)
  coef_2b = su_model2b$coefficients
  fitted_pb = model2b$fitted.values
  # c -- P(R|T)
  model2c = glm(R ~ T + X1 + D, data = cad_pseudo, family = "binomial", weight = weight_k)
  su_model2c = summary(model2c)
  coef_2c = su_model2c$coefficients
  fitted_pc = model2c$fitted.values
  # all
  coef_2 = append(coef_2, list(list(coef_2a, coef_2b, coef_2c)))
  # weight
  fitted_ps = list(fitted_pa, fitted_pb, fitted_pc)
  ys = list(model2a$y, model2b$y, model2c$y)
  p0 = (fitted_ps[[1]]^ys[[1]]) * ((1 - fitted_ps[[1]])^(1 - ys[[1]]))  # P(D)
  p1 = (fitted_ps[[2]]^ys[[2]]) * ((1 - fitted_ps[[2]])^(1 - ys[[2]]))  # P(T|D)
  p2 = (fitted_ps[[3]]^ys[[3]]) * ((1 - fitted_ps[[3]])^(1 - ys[[3]]))  # P(R|T)
  pk = p0 * p1 * p2  # P(R,T,D|X)
  weight_k[index_2] = pk[index_2] / (pk[index_2] + pk[index_3])
  weight_k[index_3] = 1 - weight_k[index_2]
  # next iteration
  t = t + 1
}
coef_2[[max_t]]
# Marginal estimates
cad_1 = cad_0 = cad
cad_1$D = 1  # set all D to 1
cad_0$D = 0  # set all D to 0
sn4x = mean(prediction(model2b, data=cad_1)$fitted * prediction(model2a, data=cad_1)$fitted) / 
  mean(prediction(model2a, data=cad_1)$fitted); sn4x      # P(T=1|D=1,X1)
sp4x = mean((1 - prediction(model2b, data=cad_0)$fitted) * (1 - prediction(model2a, data=cad_0)$fitted)) /
  mean(1 - prediction(model2a, data=cad_0)$fitted); sp4x  # P(T=0|D=0,X1)

# Compare
snsp = c(sn0,sp0,sn1,sp1,sn2,sp2,sn2x,sp2x,sn3,sp3,sn3x,sp3x,sn3a,sp3a,sn3ax,sp3ax,
         sn4,sp4,sn4x,sp4x)
snsp = matrix(snsp, ncol = 2, nrow = length(snsp)/2, byrow = T)
rownames(snsp) = c("Complete-case Analysis","B&G Count",
                   "B&G Regression","B&G Regression with Covariate",
                   "Multiple Imputation (LogReg)","Multiple Imputation with Covariate (LogReg)",
                   "Multiple Imputation (PMM)","Multiple Imputation with Covariate (PMM)",
                   "EM Algorithm MNAR","EM Algorithm MNAR with Covariate")
colnames(snsp) = c("Sensitivity","Specificity")
tbl_compare = round(snsp, 3); tbl_compare
knitr::kable(tbl_compare, "rst")
