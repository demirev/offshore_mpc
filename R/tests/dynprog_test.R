source("R/functions/dynprog.R")
source("R/functions/utils.R")
source("R/classes/shocks_util_prod.R")

# test1 <- pf_iter() # default settings - no shocks
# 
# test2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 200
# )
# test2_2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.1), 
#   psi = LogNormalShock$new(sigma = 0.1), 
#   ndraw = 200
# )
# 
# test3 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# test3_2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 1500
# )
# test3_3 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.1), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# 
# test4 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# test4_2 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 1500
# )
# test4_3 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.1), 
#   ndraw = 500
# )

# test5_1 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.5),
#   psi = LogNormalShock$new(sigma = 0.5),
#   ndraw = 200
# )
# test5_2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.5),
#   #psi = LogNormalShock$new(sigma = 0.5),
#   ndraw = 200
# )
# test5_3 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.5),
#   psi = LogNormalShock$new(sigma = 0.5),
#   ndraw = 200
# )

# test6_1 <- pf_iter(
#   psi = LogNormalShock$new(sigma = 0.0025),
#   ndraw = 400
# )
# test6_2 <- pf_iter(
#   xi = EmploymentShock$new(sigma = 0.04),
#   ndraw = 400
# )
# test6_3 <- pf_iter(
#   xi = EmploymentShock$new(sigma = 0.04),
#   psi = LogNormalShock$new(sigma = 0.025),
#   ndraw = 60
# ) 
# 
# test7_1 <- pf_iter(
#   xi = EmploymentShock$new(sigma = 0.04),
#   psi = LogNormalShock$new(sigma = 0.025),
#   ndraw = 60,
#   m_seq = discretize_m(
#     max_m = 35,
#     num_out = 130
#   ),
#   k_seq = c(26, 38, 46, 54)
# )


system.time(
test8_1 <- pf_iter(
  xi = EmploymentShock$new(sigma = 0.04),
  psi = LogNormalShock$new(sigma = 0.025),
  ndraw = 60,
  m_seq = discretize_m(
    max_m = 35,
    num_out = 130
  ),
  k_seq = 40, 
  tol = 0.025
)
)

# test8_2 <- pf_iter(
#   xi = EmploymentShock$new(sigma = 0.04),
#   psi = LogNormalShock$new(sigma = 0.025),
#   ndraw = 60,
#   m_seq = discretize_m(
#     max_m = 35,
#     num_out = 130
#   ),
#   k_seq = 40, 
#   action = test8_1$QTable$action, fit_policy = fit_spline2
# )
# 
# 
# fit_spline3 <- function(V, df_m = NULL, df_k = NULL) {
#   
#   if (is.null(df_k)) df_k <- length(unique(V$k))/3
#   if (is.null(df_m)) df_m <- length(unique(V$m))/3
#   
#   if (length(unique(V$k)) > 1) {
#     fit <- gam(
#       action ~ lo(m, span = 10) + bs(k, span = 10), 
#       data = V
#     )
#   } else {
#     fit <- gam(
#       action ~  lo(m, span = 10), 
#       data = V
#     )
#   }
#   
#   pred_func <- function(m,k){
#     if (length(m) > 1 & length(k) == 1) {
#       k <- rep(k, length(m))
#     }
#     predict(fit, data.frame(m, k, value = 0), type = "response")
#   }
#   return(pred_func)
# }
# 
# fit <- fit_spline3(test8_1$QTable, df_m = NULL)
# 
# policy_plot(fit, test8_1$QTable)

