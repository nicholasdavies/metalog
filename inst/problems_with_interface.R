library(metalog)

mlog(mu = c(0, 1))
# Error in mlog(mu = c(0, 1)) : argument "x" is missing, with no default

mlog(x = NULL, mu = c(0, 1))
# Error in metalog:::xtabulate(x) :
# Not compatible with requested type: [type=NULL; target=double].


# Works but not easy
mlog(x = NULL, y = "moments", mu = c(0, 1))
# 2-term metalog distribution with bounds (-Inf, Inf) fit using exact search by root finding
# Coefficients: 0 0.5513289
# Contains elements: a b cache method

# Unexpected behaviour - probably need to ditch y as confidence interval interpretation
mlog(x = 0.5, y = 0.75, mu = c(0, 1))
# Error in mlog(x = 0.5, y = 0.75, mu = c(0, 1)) :
#   cannot interpret y as confidence interval unless x has 2 or 3 elements
