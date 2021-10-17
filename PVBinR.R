# Info ====
# PVB Correction Tutorial
# Author: Wan Nor Arifin


# Load PVBcorrect package ====

# Install the package, if this is not yet installed, uncomment the lines below
# install.packages("devtools")  # uncomment
# devtools::install_github("wnarifin/PVBcorrect")  # uncomment
library(PVBcorrect)


# Load and explore data ====

# Data, Kosinski & Barnhart (2003)
# T = SPECT thallium test: 1 = positive, 0 = negative.
# D = CAD (disease) status: 1 = yes, 0 = no; NA = not verified.
# X1 = Gender: 1 = male, 0 = female.
# X2 = Stress mode: 1 = dipyridamole, 0 = exercise.
# X3 = Age: 1 = age ≥ 60 years, 0 = age < 60 years.

# Load data
?cad_pvb  # data info, built in data in the package

# Explore the data
str(cad_pvb)
head(cad_pvb); tail(cad_pvb)
summary(cad_pvb)

# View table
view_table(data = cad_pvb, test = "T", disease = "D")
view_table(data = cad_pvb, test = "T", disease = "D", show_unverified = TRUE)


# Obtain the accuracy estimates ====

## CCA ====

# no covariate
cca_out = acc_cca(data = cad_pvb, test = "T", disease = "D", ci = TRUE)
cca_est = cca_out$acc_results
cca_est

## EBG ====

# no covariate
ebg_out = acc_ebg(data = cad_pvb, test = "T", disease = "D", ci = TRUE, ci_type = "bca",
                  seednum = 12345, R = 999)
ebg_est = ebg_out$acc_results
ebg_est

# with covariate
ebgx_out = acc_ebg(data = cad_pvb, test = "T", disease = "D", covariate = "X3", saturated_model = TRUE,
                   ci = TRUE, ci_type = "bca", seednum = 12345, R = 999)
ebgx_est = ebgx_out$acc_results
ebgx_est

## MI ====

# no covariate
mi_out = acc_mi(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345, m = 85)
mi_est = mi_out$acc_results
mi_est

# with covariate
mix_out = acc_mi(data = cad_pvb, test = "T", disease = "D", covariate = "X3", ci = TRUE, seednum = 12345, m = 85)
mix_est = mix_out$acc_results
mix_est

## EM ====

# no covariate
# save to an R object for detailed analysis later
# for sample run, test with low R boot number, will take very long time to finish for large R
start_time = proc.time()
em_out = acc_em(data = cad_pvb, test = "T", disease = "D", ci = TRUE, ci_type = "bca",
                seednum = 12345, R = 999,
                t_max = 5000, cutoff = 0.0002)
elapsed_time = proc.time() - start_time; elapsed_time["elapsed"]  # time taken in seconds
em_est = em_out$acc_results
em_est

# with covariate, will take some time
# check the time taken to finish
startx_time = proc.time()
emx_out = acc_em(data = cad_pvb, test = "T", disease = "D", covariate = "X3", ci = TRUE, ci_type = "bca",
                 seednum = 12345, R = 999,
                 t_max = 50000, cutoff = 0.0002)  # with covariate, better set larger t_max
elapsedx_time = proc.time() - startx_time; elapsedx_time["elapsed"]  # time taken in seconds
emx_est = emx_out$acc_results
emx_est


# Print a combined table for all results ====
tbl_combined = tibble::tibble(
  Estimates = c("Sensitivity", "SE", "2.5%", "97.5%", "Specificity", "SE", "2.5%", "97.5%", 
                "PPV", "SE", "2.5%", "97.5%", "NPV", "SE", "2.5%", "97.5%"),
  CCA = c(t(cca_est["Sn", ]), t(cca_est["Sp", ]), t(cca_est["PPV", ]), t(cca_est["NPV", ])),
  EBG = c(t(ebg_est["Sn", ]), t(ebg_est["Sp", ]), t(ebg_est["PPV", ]), t(ebg_est["NPV", ])),
  EBGX = c(t(ebgx_est["Sn", ]), t(ebgx_est["Sp", ]), t(ebgx_est["PPV", ]), t(ebgx_est["NPV", ])),
  MI = c(t(mi_est["Sn", ]), t(mi_est["Sp", ]), t(mi_est["PPV", ]), t(mi_est["NPV", ])),
  MIX = c(t(mix_est["Sn", ]), t(mix_est["Sp", ]), t(mix_est["PPV", ]), t(mix_est["NPV", ])),
  EM = c(t(em_est["Sn", ]), t(em_est["Sp", ]), t(em_est["PPV", ]), t(em_est["NPV", ])),
  EMX = c(t(emx_est["Sn", ]), t(emx_est["Sp", ]), t(emx_est["PPV", ]), t(emx_est["NPV", ]))
); tbl_combined

tbl_combined_rnd = tbl_combined
tbl_combined_rnd[-1] = round(tbl_combined[-1], 3)
tbl_view = knitr::kable(tbl_combined_rnd, format = "simple",
             caption = "Comparison of accuracy estimates by PVB correction methods.")
  # format = "simple", html", "latex", "pipe", "rst"
tbl_view

# Check EM convergence ====
# all t should be < t_max for each, will return TRUE
max(em_out$boot_data$t[, 5]) < 5000  # w/out covariate, all converged
max(emx_out$boot_data$t[, 5]) < 50000  # with covariate, all converged
