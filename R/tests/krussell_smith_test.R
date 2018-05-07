# kStart = 52.95991 # 38.9
# nAgent = 10000
# testKs <- KS_Economy$new(
#   Agents = list(
#     AgentType$new(
#       M = kStart, 
#       p = 1, 
#       l = 1/0.9,
#       L = 1,
#       psi = LogNormalShock$new(sigma = 0.025), 
#       xi = EmploymentShock$new(sigma = 0.04), 
#       utilf = Iso_Elastic$new(rho = 1),
#       pol = function(m, k) return(m*0.5),
#       QTable = NULL,
#       D = 0.00625,
#       beta = 0.985,
#       n = nAgent
#     )
#   ),
#   K = kStart * (nAgent/0.9), 
#   P = 1, 
#   FF = Cobb_Douglas$new(alpha = 0.36),
#   PSI = NoShock$new(mu = 1),
#   XI  = NoShock$new(mu = 1),
#   delta = 0.025 
# )
# 
# #testKs$stohastic_optimization(k = kStart) #52.95991
# testKs$update_policies(k = testKs$K / testKs$L)
# testKs$wealth_ss(verbose = T)
