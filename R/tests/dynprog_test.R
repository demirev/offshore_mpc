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