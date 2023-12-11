library(sox)

# overlapping groups - graph proximal
# grouping structure
# ## group
# grp <- matrix(c(0, 0, 0, 0, 0,
#                 0, 0, 0, 0, 0,
#                 1, 1, 0, 0, 0,
#                 0, 0, 0, 0, 0,
#                 0, 1, 0, 1, 0),
#               ncol = 5, byrow = TRUE)
# ## group_variable
# grp.var <- matrix(c(1, 0, 0, 0, 0, #A1
#                     1, 0, 0, 0, 0, #A2
#                     0, 0, 0, 1, 0, #C1
#                     0, 0, 0, 1, 0, #C2
#                     0, 1, 0, 0, 0, #B
#                     0, 0, 1, 0, 0, #A1B
#                     0, 0, 1, 0, 0, #A2B
#                     0, 0, 0, 0, 1, #C1B
#                     0, 0, 0, 0, 1  #C2B
# ), ncol = 5, byrow = TRUE)

x <- as.matrix(sim[, c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B")])
lam.seq <- exp(seq(log(1e0), log(1e-3), length.out = 20))

# Variables:
## 1: A1
## 2: A2
## 3: C1
## 4: C2
## 5: B
## 6: A1B
## 7: A2B
## 8: C1B
## 9: C2B

# Overlapping groups:
## g1: A1, A2, A1B, A2B
## g2: B, A1B, A2B, C1B, C2B
## g3: A1B, A2B
## g4: C1, C2, C1B, C2B
## g5: C1B, C2B

# So..
overlapping.groups <- list(c(1, 2, 6, 7),
                           c(5, 6, 7, 8, 9),
                           c(6, 7),
                           c(3, 4, 8, 9),
                           c(8, 9))
pars.overlapping <- overlap_structure(overlapping.groups)

# fit
fit.overlapping <- sox(
  x = x,
  ID = sim$Id,
  time = sim$Start,
  time2 = sim$Stop,
  event = sim$Event,
  penalty = "overlapping",
  lambda = lam.seq,
  group = pars.overlapping$groups,
  group_variable = pars.overlapping$groups_var,
  penalty_weights = pars.overlapping$group_weights,
  tol = 1e-4,
  maxit = 1e3,
  verbose = FALSE
)
plot(fit.overlapping)

# cv
cv.overlapping <- sox_cv(
  x = x,
  ID = sim$Id,
  time = sim$Start,
  time2 = sim$Stop,
  event = sim$Event,
  penalty = "overlapping",
  lambda = lam.seq,
  group = pars.overlapping$groups,
  group_variable = pars.overlapping$groups_var,
  penalty_weights = pars.overlapping$group_weights,
  nfolds = 5,
  tol = 1e-4,
  maxit = 1e3,
  verbose = FALSE
)
plot(cv.overlapping$sox.fit)
plot(cv.overlapping)
plot(cv.overlapping, type = "solution-path")




# Variables:
## 1: A1
## 2: A2
## 3: C1
## 4: C2
## 5: B
## 6: A1B
## 7: A2B
## 8: C1B
## 9: C2B

# Nested groups (misspecified, for the demonstration of the software only.)
## g1: A1, A2, C1, C2, B, A1B, A2B, C1B, C2B
## g2: A1B, A2B, A1B, A2B
## g3: C1, C2, C1B, C2B
## g4: 1
## g5: 2
## ...
## G12: 9

# So...
nested.groups <- list(1:9,
                      c(1, 2, 6, 7),
                      c(3, 4, 8, 9),
                      1, 2, 3, 4, 5, 6, 7, 8, 9)

pars.nested <- nested_structure(nested.groups)

fit.nested <- sox(
  x = x,
  ID = sim$Id,
  time = sim$Start,
  time2 = sim$Stop,
  event = sim$Event,
  penalty = "nested",
  lambda = lam.seq,
  group = pars.nested$groups,
  own_variable = pars.nested$own_variables,
  no_own_variable = pars.nested$N_own_variables,
  penalty_weights = pars.nested$group_weights,
  tol = 1e-4,
  maxit = 1e3,
  verbose = FALSE
)
plot(fit.nested)

# cv
cv.nested <- sox_cv(
  x = x,
  ID = sim$Id,
  time = sim$Start,
  time2 = sim$Stop,
  event = sim$Event,
  penalty = "nested",
  lambda = lam.seq,
  group = pars.nested$groups,
  own_variable = pars.nested$own_variables,
  no_own_variable = pars.nested$N_own_variables,
  penalty_weights = pars.nested$group_weights,
  nfolds = 5,
  tol = 1e-4,
  maxit = 1e3,
  verbose = FALSE
)
plot(cv.nested$sox.fit)
plot(cv.nested)
plot(cv.nested, type = "solution-path")
