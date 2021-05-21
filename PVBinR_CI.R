# PVB Correction -- with 95% CI

# Data, 2003a Kosinski & Barnhart
# X1 = Gender: 1 = male, 0 = female.
# X2 = Stress mode: 1 = dipyridamole, 0 = exercise.
# X3 = Age: 1 = age ≥ 60 years, 0 = age < 60 years.
# T = SPECT thallium test: 1 = positive, 0 = negative.
# D = CAD (disease) status: 1 = yes, 0 = no; NA = not verified.
# May try with other data, modify variable names accordingly.
data = read.csv("2003kosinski_cad.csv")
str(data)

# Verification status
data$V = 1  # V = Verified: 1 = yes, 0 = no
data[is.na(data$D), "V"] = 0

# Percent missing
mean(is.na(data$D))*100

# Library
library(mice)
library(boot)
library(prediction)

# Session Info
sessionInfo()

# Functions
# For use with CCA, MI
# Sn with variance
sn_calc_with_var = function(data) {
  d = data
  tbl = table(d$T, d$D)
  sn = tbl[2,2]/sum(tbl[,2])
  var = tbl[2,2]*tbl[1,2] / sum(tbl[,2])^3
  return(cbind(sn = sn, sn_var = var))
}
# Sp with variance
sp_calc_with_var = function(data) {
  d = data
  tbl = table(d$T, d$D)
  sp = tbl[1,1]/sum(tbl[,1])
  var = tbl[1,1]*tbl[2,1] / sum(tbl[,1])^3
  return(cbind(sp = sp, sp_var = var))
}
# For use with bootstrap
# Sn Sp for EBG using GLM
snsp_ebg = function(data, indices) {
  d = data[indices,]
  # data = complete case only for modeling
  fit = glm(D ~ T, data = d, family = "binomial")
  preds = prediction(fit, data_minus)$fitted  # Predicted for complete & incomplete data
  sn = sum(data_minus$T*preds) / sum(preds)  # P(T=1|D=1)
  sp = sum((1-data_minus$T)*(1-preds)) / sum(1-preds)  # P(T=0|D=0)
  return(cbind(sn, sp))
}
# Sn Sp for EBG using GLM, one covariate
snsp_ebgx = function(data, indices) {
  d = data[indices,]
  # data = complete case only for modeling
  # fit = glm(D ~ T + X1, data = d, family = "binomial")  # Alonzo, unsaturated
  fit = glm(D ~ T*X1, data = d, family = "binomial")  # saturated
  preds = prediction(fit, data_minus)$fitted  # Predicted for complete & incomplete data
  sn = sum(data_minus$T*preds) / sum(preds)  # P(T=1|D=1,X)
  sp = sum((1-data_minus$T)*(1-preds)) / sum(1-preds)  # P(T=0|D=0,X)
  return(cbind(sn, sp))
}

# Other options
seednum = 12345  # seednum for bootstrap and MI
m = 85  # number of imputed data sets for MI
b = 999  # number of bootstrap

# Complete-Case Analysis
# sn
sn_cc_all = sn_calc_with_var(na.omit(data))
sn_cc = sn_cc_all[,"sn"]
sn_cc_se = sqrt(sn_cc_all[,"sn_var"])
sn_cc_ci = sn_cc + c(-1,1) * qnorm(1 - 0.05/2) * sn_cc_se
c(sn_cc, sn_cc_se, sn_cc_ci)
# sp
sp_cc_all = sp_calc_with_var(na.omit(data))
sp_cc = sp_cc_all[,"sp"]
sp_cc_se = sqrt(sp_cc_all[,"sp_var"])
sp_cc_ci = sp_cc + c(-1,1) * qnorm(1 - 0.05/2) * sp_cc_se
c(sp_cc, sp_cc_se, sp_cc_ci)

# B&G, Count Method in Begg & Greenes 1983
# Not implemented here, use regression method instead, easier to setup with bootstrap

# Recode NA to -1 for ease of calculation for B&G
data_minus = data
data_minus[is.na(data_minus$D), "D"] = -1  # recode NA as -1, this will be -1*0=0 for unverified. Else NA*0=NA.

# B&G, Regression Method in Alonzo 2005
set.seed(seednum)
snsp_ebg_boot_data = boot(data = na.omit(data), statistic = snsp_ebg, R = b)
snsp_ebg_boot = snsp_ebg_boot_data$t0
snsp_ebg_boot_se = apply(snsp_ebg_boot_data$t, 2, sd)
# sn
sn_ebg_boot = snsp_ebg_boot[1]
sn_ebg_boot_se = snsp_ebg_boot_se[1]
sn_ebg_boot_ci_data = boot.ci(snsp_ebg_boot_data, type = "bca", index = 1)
sn_ebg_boot_ci = sn_ebg_boot_ci_data$bca[4:5]
c(sn_ebg_boot, sn_ebg_boot_se, sn_ebg_boot_ci)
# sp
sp_ebg_boot = snsp_ebg_boot[2]
sp_ebg_boot_se = snsp_ebg_boot_se[2]
sp_ebg_boot_ci_data = boot.ci(snsp_ebg_boot_data, type = "bca", index = 2)
sp_ebg_boot_ci = sp_ebg_boot_ci_data$bca[4:5]
c(sp_ebg_boot, sp_ebg_boot_se, sp_ebg_boot_ci)

# B&G, Regression Method in Alonzo 2005, with one covariate
set.seed(seednum)
snsp_ebgx_boot_data = boot(data = na.omit(data), statistic = snsp_ebgx, R = b)
snsp_ebgx_boot = snsp_ebgx_boot_data$t0
snsp_ebgx_boot_se = apply(snsp_ebgx_boot_data$t, 2, sd)
# sn
sn_ebgx_boot = snsp_ebgx_boot[1]
sn_ebgx_boot_se = snsp_ebgx_boot_se[1]
sn_ebgx_boot_ci_data = boot.ci(snsp_ebgx_boot_data, type = "bca", index = 1)
sn_ebgx_boot_ci = sn_ebgx_boot_ci_data$bca[4:5]
c(sn_ebgx_boot, sn_ebgx_boot_se, sn_ebgx_boot_ci)
# sp
sp_ebgx_boot = snsp_ebgx_boot[2]
sp_ebgx_boot_se = snsp_ebgx_boot_se[2]
sp_ebgx_boot_ci_data = boot.ci(snsp_ebgx_boot_data, type = "bca", index = 2)
sp_ebgx_boot_ci = sp_ebgx_boot_ci_data$bca[4:5]
c(sp_ebgx_boot, sp_ebgx_boot_se, sp_ebgx_boot_ci)

# Multiple Imputation, Logistic Regression

# Setup - no covariate
data1 = data[, c("T","D")]  # select only T & D, else MI use all variables for imputation
data1$D = as.factor(data1$D)  # D as factor

data_mi_logreg = mice(data1, m = m, method = "logreg", seed = seednum, print = F)  # MIDS class
data_mi_logreg_data = complete(data_mi_logreg, "all")  # imputed data

# Imputed Sn, Sp with var each
sn_mi_logreg_imp = t(sapply(data_mi_logreg_data, sn_calc_with_var))
colnames(sn_mi_logreg_imp) = c("sn", "sn_var")
sp_mi_logreg_imp = t(sapply(data_mi_logreg_data, sp_calc_with_var))
colnames(sp_mi_logreg_imp) = c("sp", "sp_var")

# Pooled Sn, Sp
sn_mi_logreg_pool = pool.scalar(Q = sn_mi_logreg_imp[,"sn"], U = sn_mi_logreg_imp[,"sn_var"])
sp_mi_logreg_pool = pool.scalar(Q = sp_mi_logreg_imp[,"sp"], U = sp_mi_logreg_imp[,"sp_var"])

# Point & CI
# sn
sn_mi_logreg = sn_mi_logreg_pool$qbar
# z-dist
# sn_mi_logreg_ci = sn_mi_logreg_pool$qbar + c(-1,1) * qnorm(1 - 0.05/2) * sqrt(sn_mi_logreg_pool$t)
# t-dist, Rubin's rule
v_sn = (m-1)*(1+(sn_mi_logreg_pool$ubar/((1+1/m)*sn_mi_logreg_pool$b)))^2
t_sn = qt((1-0.05/2), v_sn)
sn_mi_logreg_se = sqrt(sn_mi_logreg_pool$t)
sn_mi_logreg_ci = sn_mi_logreg_pool$qbar + c(-1,1) * t_sn * sn_mi_logreg_se
c(sn_mi_logreg, sn_mi_logreg_se, sn_mi_logreg_ci)
# sp
sp_mi_logreg = sp_mi_logreg_pool$qbar
# z-distr
# sp_mi_logreg_ci = sp_mi_logreg_pool$qbar + c(-1,1) * qnorm(1 - 0.05/2) * sqrt(sp_mi_logreg_pool$t)
# t-dist, Rubin's rule
v_sp = (m-1)*(1+(sp_mi_logreg_pool$ubar/((1+1/m)*sp_mi_logreg_pool$b)))^2
t_sp = qt((1-0.05/2), v_sp)
sp_mi_logreg_se = sqrt(sp_mi_logreg_pool$t)
sp_mi_logreg_ci = sp_mi_logreg_pool$qbar + c(-1,1) * t_sp * sp_mi_logreg_se
c(sp_mi_logreg, sp_mi_logreg_se, sp_mi_logreg_ci)

# Multiple Imputation, Logistic Regression, with one covariate

# Setup - no covariate
data1 = data[, c("T","D","X1")]  # select only T, D & X1, else MI use all variables for imputation
data1$D = as.factor(data1$D)  # D as factor

data_mi_logregx = mice(data1, m = m, method = "logreg", seed = seednum, print = F)  # MIDS class
data_mi_logregx_data = complete(data_mi_logregx, "all")  # imputed data

# Imputed Sn, Sp with var each
sn_mi_logregx_imp = t(sapply(data_mi_logregx_data, sn_calc_with_var))
colnames(sn_mi_logregx_imp) = c("sn", "sn_var")
sp_mi_logregx_imp = t(sapply(data_mi_logregx_data, sp_calc_with_var))
colnames(sp_mi_logregx_imp) = c("sp", "sp_var")

# Pooled Sn, Sp
sn_mi_logregx_pool = pool.scalar(Q = sn_mi_logregx_imp[,"sn"], U = sn_mi_logregx_imp[,"sn_var"])
sp_mi_logregx_pool = pool.scalar(Q = sp_mi_logregx_imp[,"sp"], U = sp_mi_logregx_imp[,"sp_var"])

# Point & CI
# sn
sn_mi_logregx = sn_mi_logregx_pool$qbar
# z-dist
# sn_mi_logregx_ci = sn_mi_logregx_pool$qbar + c(-1,1) * qnorm(1 - 0.05/2) * sqrt(sn_mi_logregx_pool$t)
# t-dist, Rubin's rule
v_sn = (m-1)*(1+(sn_mi_logregx_pool$ubar/((1+1/m)*sn_mi_logregx_pool$b)))^2
t_sn = qt((1-0.05/2), v_sn)
sn_mi_logregx_se = sqrt(sn_mi_logregx_pool$t)
sn_mi_logregx_ci = sn_mi_logregx_pool$qbar + c(-1,1) * t_sn * sn_mi_logregx_se
c(sn_mi_logregx, sn_mi_logregx_se, sn_mi_logregx_ci)
# sp
sp_mi_logregx = sp_mi_logregx_pool$qbar
# z-distr
# sp_mi_logregx_ci = sp_mi_logregx_pool$qbar + c(-1,1) * qnorm(1 - 0.05/2) * sqrt(sp_mi_logregx_pool$t)
# t-dist, Rubin's rule
v_sp = (m-1)*(1+(sp_mi_logregx_pool$ubar/((1+1/m)*sp_mi_logregx_pool$b)))^2
t_sp = qt((1-0.05/2), v_sp)
sp_mi_logregx_se = sqrt(sp_mi_logregx_pool$t)
sp_mi_logregx_ci = sp_mi_logregx_pool$qbar + c(-1,1) * t_sp * sp_mi_logregx_se
c(sp_mi_logregx, sp_mi_logregx_se, sp_mi_logregx_ci)

# EM Algorithm, MNAR
# Kosinski & Barnhart, 2003
# Components: a = disease, b = diagnostic, c = missing data mechanism

# EM Algorithm Function
pvb_em = function(data_pseudo, t_max, cutoff, a, b, c, index_1, index_2, index_3) {
  # # Initialize values
  coef = vector("list", 0)  # empty vector list of coef

  # EM Iteration
  t = 1
  while (t < t_max + 1) {
    # a -- P(D|X)
    model_a = glm(a, data = data_pseudo, family = "binomial", weight = weight_k)
    su_model_a = summary(model_a)
    coef_a = su_model_a$coefficients
    fitted_pa = model_a$fitted.values
    # b -- P(T|D,X)
    model_b = glm(b, data = data_pseudo, family = "binomial", weight = weight_k)
    su_model_b = summary(model_b)
    coef_b = su_model_b$coefficients
    fitted_pb = model_b$fitted.values
    # c -- P(V|T,X,D)
    model_c = glm(c, data = data_pseudo, family = "binomial", weight = weight_k)
    su_model_c = summary(model_c)
    coef_c = su_model_c$coefficients
    fitted_pc = model_c$fitted.values
    # all
    coef = append(coef, list(list(coef_a, coef_b, coef_c)))
    # weight
    fitted_ps = list(fitted_pa, fitted_pb, fitted_pc)
    ys = list(model_a$y, model_b$y, model_c$y)
    p0 = (fitted_ps[[1]]^ys[[1]]) * ((1 - fitted_ps[[1]])^(1 - ys[[1]]))  # P(D|X)
    p1 = (fitted_ps[[2]]^ys[[2]]) * ((1 - fitted_ps[[2]])^(1 - ys[[2]]))  # P(T|D,X)
    p2 = (fitted_ps[[3]]^ys[[3]]) * ((1 - fitted_ps[[3]])^(1 - ys[[3]]))  # P(V|T,X,D)
    pk = p0 * p1 * p2  # P(V,T,D)
    weight_k[index_2] <<- pk[index_2] / (pk[index_2] + pk[index_3])
    weight_k[index_3] <<- 1 - weight_k[index_2]
    # check if change in coef < 0.001 in abs value
    if (t < 2) {
      diffs_t = cutoff + 1
    } else {  # diffs in coef
      nrow_a = nrow(coef[[t]][[1]])
      nrow_b = nrow(coef[[t]][[2]])
      nrow_c = nrow(coef[[t]][[3]])
      diffs_t_a = coef[[t]][[1]][1:nrow_a] - coef[[t-1]][[1]][1:nrow_a]  # comp a
      diffs_t_b = coef[[t]][[2]][1:nrow_b] - coef[[t-1]][[2]][1:nrow_b]  # comp b
      diffs_t_c = coef[[t]][[3]][1:nrow_c] - coef[[t-1]][[3]][1:nrow_c]  # comp c
      diffs_t = c(diffs_t_a, diffs_t_b, diffs_t_c)
    }
    # early stopping rule
    if (all(abs(diffs_t) < cutoff)) {
      cat(paste("Iteration stops early at t =", t, "because all changes <", cutoff, "\n"))
      break
    }
    # if early stopping rule is not met, run until t max
    else {
      # printing frequency
      if(t %% 50 == 0) {cat(paste("Current t =", t, "\n"))}
      t = t + 1
    }
  }
  
  # Outputs
  out = list(model_a = model_a, model_b = model_b, model_c = model_c,
             diffs_t = diffs_t, t = t-1, coef = coef)
}

# Sn Sp Function for EM
# basic, no marginal, faster
snsp_em0 = function(data_verified, indices, t_max, cutoff, data_unverified, a, b, c, index_1, index_2, index_3) {
  d = data_verified[indices, ]
  data_pseudo = rbind(d, data_unverified)
  # the EM run
  pvb_em_out = pvb_em(data_pseudo, t_max, cutoff, a, b, c, index_1, index_2, index_3)
  sn = predict(pvb_em_out$model_b, list(D=1), type="response")  # P(T=1|D=1)
  sp = 1 - predict(pvb_em_out$model_b, list(D=0), type="response")  # P(T=0|D=0)
  cat(paste("Boot Iteration =", counter, "\n"))
  cat("==================")
  cat("\n")
  counter <<- counter + 1  # <<- so that counter won't reset with each boot iter
  return(cbind(sn, sp, t = pvb_em_out$t-1))
}

# Sn Sp Function for EM
# with marginal
snsp_em = function(data_verified, indices, t_max, cutoff, data_unverified, a, b, c, index_1, index_2, index_3) {
  d = data_verified[indices, ]
  data_pseudo = rbind(d, data_unverified)
  # the EM run
  pvb_em_out = pvb_em(data_pseudo, t_max, cutoff, a, b, c, index_1, index_2, index_3)
  # prepare data for marginal estimate
  data_1 = data_0 = data_pseudo[c(index_1, index_2),]
  data_1$D = 1; data_0$D = 0
  # calculate marginal estimates
  sn = sum(prediction(pvb_em_out$model_b, data=data_1)$fitted * prediction(pvb_em_out$model_a, data=data_1)$fitted) / 
    sum(prediction(pvb_em_out$model_a, data=data_1)$fitted)  # P(T=1|D=1,X)
  sp = sum((1 - prediction(pvb_em_out$model_b, data=data_0)$fitted) * (1 - prediction(pvb_em_out$model_a, data=data_0)$fitted)) /
    sum(1 - prediction(pvb_em_out$model_a, data=data_0)$fitted)  # P(T=0|D=0,X)
  cat(paste("Boot Iteration =", counter, "\n"))
  cat("==================")
  cat("\n")
  counter <<- counter + 1  # <<- so that counter won't reset with each boot iter
  # print(indices)  only enable if used w/in boot function, to display boot samples
  return(cbind(sn, sp, t = pvb_em_out$t))
}

# Prepare Pseudo-data for EM boot
data_verified = data[data$V == 1, ]
data_unverified = data[data$V == 0, ]  # 1U, unverified U observation
data_unverified = rbind(data_unverified, data[data$V == 0, ])  # 2U, replicate U observation
# create 0, 1 for 2U rows
data_unverified[1:nrow(data[data$V == 0, ]), "D"] = 0        # 1st U rows
data_unverified[(nrow(data[data$V == 0, ])+1):nrow(data_unverified), "D"] = 1  # 2nd U rows
# index k
index_1 = which(data$V == 1)  # verified
index_2 = which(data$V == 0)  # unverified U
index_3 = max(index_2)+(1:length(index_2))  # unverified 2U

# MNAR: a = D~1, b = T~D, c = V~T+D
# bootstrap
# Setup
t_max = 500  # max iteration, t
cutoff = 0.0002  # early stopping rule if changes in coef < cutoff
ptm0 = proc.time()
set.seed(seednum)
counter = 0  # start counter at 0
weight_k = rep(1, max(index_3))  # initialize weight = 1 for all
snsp_em_boot_data = boot(data = data_verified, statistic = snsp_em0, R = b, t_max = t_max, cutoff = cutoff,
                         data_unverified = data_unverified, a = "D~1", b = "T~D", c = "V~T+D",
                         index_1 = index_1, index_2 = index_2, index_3 = index_3)
ptm1 = proc.time() - ptm0; ptm1  # view time taken in seconds
# bootstrap gets extra +1 bcs first boot run is original data
snsp_em_boot = snsp_em_boot_data$t0
snsp_em_boot_se = apply(snsp_em_boot_data$t, 2, sd)
# sn
sn_em_boot = snsp_em_boot[1]
sn_em_boot_se = snsp_em_boot_se[1]
sn_em_boot_ci_data = boot.ci(snsp_em_boot_data, type = "bca", index = 1)
sn_em_boot_ci = sn_em_boot_ci_data$bca[4:5]
# sp
sp_em_boot = snsp_em_boot[2]
sp_em_boot_se = snsp_em_boot_se[2]
sp_em_boot_ci_data = boot.ci(snsp_em_boot_data, type = "bca", index = 2)
sp_em_boot_ci = sp_em_boot_ci_data$bca[4:5]

# MNAR with Covariate: a = D~X1, b = T~D+X1, c = V~T+X1+D
# bootstrap
# Setup
t_max = 5000  # max iteration, t, needs more iter with addition of X1
cutoff = 0.0002  # early stopping rule if changes in coef < cutoff
ptm0 = proc.time()
set.seed(seednum)
counter = 0  # start counter at 0
weight_k = rep(1, max(index_3))  # initialize weight = 1 for all
snsp_emx_boot_data = boot(data = data_verified, statistic = snsp_em, R = b, t_max = t_max, cutoff = cutoff,
                         data_unverified = data_unverified, a = "D~X1", b = "T~D+X1", c = "V~T+X1+D",
                         index_1 = index_1, index_2 = index_2, index_3 = index_3)
ptm2 = proc.time() - ptm0; ptm2  # view time taken in seconds
# bootstrap gets extra +1 bcs first boot run is original data
snsp_emx_boot = snsp_emx_boot_data$t0
snsp_emx_boot_se = apply(snsp_emx_boot_data$t, 2, sd)
# sn
sn_emx_boot = snsp_emx_boot[1]
sn_emx_boot_se = snsp_emx_boot_se[1]
sn_emx_boot_ci_data = boot.ci(snsp_emx_boot_data, type = "bca", index = 1)
sn_emx_boot_ci = sn_emx_boot_ci_data$bca[4:5]
# sp
sp_emx_boot = snsp_emx_boot[2]
sp_emx_boot_se = snsp_emx_boot_se[2]
sp_emx_boot_ci_data = boot.ci(snsp_emx_boot_data, type = "bca", index = 2)
sp_emx_boot_ci = sp_emx_boot_ci_data$bca[4:5]

# Compare all methods
snsp = c(sn_cc, sn_cc_se, sn_cc_ci, sp_cc, sp_cc_se, sp_cc_ci,
         sn_ebg_boot, sn_ebg_boot_se, sn_ebg_boot_ci, sp_ebg_boot, sp_ebg_boot_se, sp_ebg_boot_ci,
         sn_ebgx_boot, sn_ebgx_boot_se, sn_ebgx_boot_ci, sp_ebgx_boot, sp_ebgx_boot_se, sp_ebg_boot_ci,
         sn_mi_logreg, sn_mi_logreg_se, sn_mi_logreg_ci, sp_mi_logreg, sp_mi_logreg_se, sp_mi_logreg_ci,
         sn_mi_logregx, sn_mi_logregx_se, sn_mi_logregx_ci, sp_mi_logregx, sp_mi_logregx_se, sp_mi_logregx_ci,
         sn_em_boot, sn_em_boot_se, sn_em_boot_ci, sp_em_boot, sp_em_boot_se, sp_em_boot_ci,
         sn_emx_boot, sn_emx_boot_se, sn_emx_boot_ci, sp_emx_boot, sp_emx_boot_se, sp_emx_boot_ci)
snsp = matrix(snsp, ncol = 8, nrow = length(snsp)/8, byrow = T)
rownames(snsp) = c("Complete-case Analysis",
                   "EBG Regression - Bootstrap CI",
                   "EBG Regression with Covariate - Bootstrap CI",
                   "Multiple Imputation (LogReg)",
                   "Multiple Imputation with Covariate (LogReg)",
                   "EM Alg. MNAR - Bootstrap CI",
                   "EM Alg. MNAR with Covariate - Bootstrap CI")
colnames(snsp) = c("Sen", "SE", "2.5%", "97.5%", "Spe", "SE", "2.5%", "97.5%")
tbl_compare = round(snsp, 3)
knitr::kable(tbl_compare, "simple")  # "simple", html", "latex", "pipe", "rst"

# Save to text files for reference
time_tag = format(Sys.time(), "%Y%m%d_%H%M%S")
write.table(tbl_compare, file = paste0("pvb_methods_comparison_x_df_", time_tag, ".txt"))
capture.output(knitr::kable(tbl_compare, "simple"), file = paste0("pvb_methods_comparison_x_formatted_", time_tag, ".txt"))
