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


# So. I think we should have params to mlog x, y, mu, bounds, n, in that order.
# Also have an additional param conf, which if specified fills in y.
# Also allow y = "data" but not y = "bin". The method for interpreting y = "data"
# should be specified by "method".
# Also allow y = "coeff". However, y = "data" and y = "coeff" cannot coexist with
# mu or n.
# Method = "auto" should have different behaviour. I think there should be a method = "approx"
# as well.

# Maybe instead of allowing y to be "data" or "coeff", scratch that and have
# mlog_coeff and mlog_fit. Check all usage of "data" and "coeff" to see if
# this is feasible.


# General methods are like "exact" or "match" -- does exact search;
# "valid" or "approx" -- attempts to enforce validity;


# Moment matching approach for fitting data?
# This doesn't work...
hist(quakes$lat) # unimodal
hist(quakes$long) # bimodal

x = quakes$lat

# Calculate moments (for mlog) - mean, sd, then kth standardized moment
mu = rep(mean(x), 9)
mu[2] = var(x)
for (i in 3:9) {
    mu[i] = e1071::moment(x, order = i, center = TRUE) / sd(x)^i
}

plot(mlog(x = NULL, y = "moments", mu = mu[1:3]))
mlog(x = NULL, y = "moments", mu = mu[1:8], method = "cvx")
mlog(x = NULL, y = "moments", mu = mu[1:8], method = "rf")

# Quantile matching?
# Works for some n but not others.
n = 6
y = (1:n - 0.5)/n
q = quantile(x, y)

# v different depending upon whether bounds are specified, n, etc
ml = mlog(q, y, bounds = c(-50, 0), method = "exact")

curve(dmlog(x, ml)*length(x)*12, from = -50, to = 10)
hist(x, add = TRUE, breaks = -50:10)
