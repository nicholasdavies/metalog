library(metalog)
library(ggplot2)
library(data.table)
library(gsl)

testfunc(1, 1)
testfunc(1, 2)
testfunc(1, 3)
testfunc(1, 4)
testfunc(1, 5)
testfunc(1, 6)
testfunc(1, 7)
testfunc(1, 8)
testfunc(1, 9)
testfunc(1, 10)
testfunc(1, 11)
testfunc(1, 12)

k=16
hyperg_2F1(1, k+2, k+3, -1, TRUE, FALSE)

plot(mlog(c(1, 2, 3), c(0.1, 0.5, 0.9), n = 3))
plot(mlog(c(1, 2, 10), c(0.1, 0.5, 0.9), n = 4, method = "LP"))

# Test of refine = TRUE and its impact on building time and dist func time
mlR = mlog(quakes$lat, "data", n = 20, refine = TRUE)
mlU = mlog(quakes$lat, "data", n = 20, refine = FALSE)

# what is kth component of mean for k even & k >= 6
a = 1
ifunc = function(y) a * (y-0.5)^(k/2-1) * log(y/(1-y));

for (k in seq(8, 32, by = 4)) {
    cat("k = ", k, ": ", 1/integrate(ifunc, lower = 0, upper = 1)$value, "\n", sep = "")
}

rmetalog::metalog(x, probs = y, term_lower_bound = 15, term_limit = 15, fit_method = "LP")

microbenchmark::microbenchmark(
    pmlog(x, mlR),
    pmlog(x, mlU),
    times = 1000,
    setup = x <- runif(1000, min(quakes$lat), max(quakes$lat))
)


dummy = pmax(0, ceiling(rnorm(2000, 20, 5)))

dt = data.table(x = dummy);
microbenchmark::microbenchmark(
    table(dummy),
    tabulate(dummy),
    dt[, .N, by = x],
    metalog:::xtabulate(dummy),
    times = 1000L
)

ml = mlog(dummy, "data", bounds = c(0, Inf), n = 5)
plot(ml)
print(ml)

x = seq(-2, 2, by = 0.05)
y = pnorm(x)
d = dnorm(x)

ml = mlog(x, y, n = 5)
plot(ml)

plot(x, y)

do_breaks = function(v)
{
    b = hist(v, plot = FALSE, right = FALSE)$breaks;
    bw = b[2] - b[1];
    if (any(v == b[length(b)])) {
        b = c(b, b[length(b)] + bw);
    }
    return (b)
}

vec = c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,6,6,6,6,6,6,6,6,6,6,8,8,8,9)
vec = rnorm(50)
hist(vec, breaks = do_breaks(vec), include.lowest = FALSE, right = FALSE)

ml = mlog(c(1,2,2.5,3,3.33,3.67,4,4.5,5), (1:9 - 0.5) / 9, n = 8)

ggplot(ml$cache) +
    geom_point(aes(M, y))

ggplot(ml$cache) +
    geom_point(aes(M, m))
    # scale_y_log10()

# Histogram check
xpts = rnorm(10000)
histy = hist(xpts, include.lowest = FALSE)

ml = mlog(xpts, "bin", n = 6)



aic = function(ml, xpts)
{
    2 * length(ml$a) - 2 * sum(dmlog(xpts, ml, log = TRUE))
}

plot(2:20, c(
    aic(mlog(histy, n = 2), xpts),
    aic(mlog(histy, n = 3), xpts),
    aic(mlog(histy, n = 4), xpts),
    aic(mlog(histy, n = 5), xpts),
    aic(mlog(histy, n = 6), xpts),
    aic(mlog(histy, n = 7), xpts),
    aic(mlog(histy, n = 8), xpts),
    aic(mlog(histy, n = 9), xpts),
    aic(mlog(histy, n = 10), xpts),
    aic(mlog(histy, n = 11), xpts),
    aic(mlog(histy, n = 12), xpts),
    aic(mlog(histy, n = 13), xpts),
    aic(mlog(histy, n = 14), xpts), # investigate this...
    aic(mlog(histy, n = 15), xpts),
    aic(mlog(histy, n = 16), xpts),
    aic(mlog(histy, n = 17), xpts),
    aic(mlog(histy, n = 18), xpts),
    aic(mlog(histy, n = 19), xpts),
    aic(mlog(histy, n = 20), xpts)
))

ml = mlog(histy, n = 14)
www = dmlog(xpts, ml, log = TRUE)
which(!is.finite(www))
xpts[317]
dmlog(xpts[317], mlog(histy, n = 14))

ggplot(ml$cache) +
    geom_point(aes(M, y))

ggplot(ml$cache) +
    geom_point(aes(y, m))

ggplot(ml$cache) +
    geom_col(data = data.frame(histy$mids, histy$density), aes(histy.mids, histy.density)) +
    geom_point(aes(M, m))



qmlog(1:10/11, ml)
qmlog(-0.000001, ml)

microbenchmark(
    pmlog(-2.5, ml),
    pnorm(-2.5),
    times = 10000
)


vars = rmlog(0, ml)
hist(vars)
vars

Q = 0
Q = runif(1000, -5, 5)

microbenchmark(
    pnorm(Q),
    pmlog(Q, ml),
    rmetalog::pmetalog(ml_old, Q, 5),
    times = 16L,
    unit = "ms"
)
bb=.Last.value
autoplot(bb)

bb2 = microbenchmark(
    dnorm(Q),
    dmlog(Q, ml),
    rmetalog::dmetalog(ml_old, Q, 5),
    times = 16L,
    unit = "ms"
)
autoplot(bb2)

bb3 = microbenchmark(
    rnorm(100000),
    rmlog(100000, ml),
    rmetalog::rmetalog(ml_old, 100000, 5),
    times = 1000L
)
bb3
autoplot(bb3)



Q = runif(1000, -5, 5)

microbenchmark(
    pmlog(runif(1000, -5, 5), ml),
    pmlog(runif(1000, -5, 5), ml),
    times = 10000L
)


x = seq(-2, 2, by = 0.05)
y = pnorm(x)
ml_old = rmetalog::metalog(x, probs = y, term_lower_bound = 5, term_limit = 5)
rmetalog::pmetalog(ml_old, -2.5, 5)

rmetalog::dmetalog(ml_old, -2.5, 5)
dmlog(-2.5, ml)

h = hist(rnorm(1000))
class(h) == "histogram"
is(h, "histogram")


# Fit bounded metalog distribution by mean and two quantiles, nloptr


library(nloptr)


eval_f = function(a, target_mean, target_x, target_y)
{
    current_mean = integrate(metalog:::metalog_M, 0, 1, a = a, n = length(a), bl = 0, bu = 1)$value;
    current_x = metalog:::metalog_M(target_y, a = a, n = length(a), bl = 0, bu = 1);
    sum( (c(target_mean, target_x) - c(current_mean, current_x))^2 )
}

eval_g = function(a, target_mean, target_x, target_y)
{
    ty = seq(0.001, 0.999, by = 0.001)
    log_odds_ty = log(ty / (1 - ty));
    tym0.5 = ty - 0.5;
    tycinv = 1 / (ty * (1 - ty));

    mat = metalog:::LP_feas_coeffs(n = length(a), m = length(ty), log_odds_ty, tym0.5, tycinv);
    -mat %*% a
}


x0 = c(0, 0, 0)
opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1.0e-18, maxeval = 10000) # NB no constraints
opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1.0e-18, maxeval = 10000) # Allows constraints


res = nloptr::nloptr(x0 = x0, eval_f = eval_f, eval_g_ineq = eval_g,
    opts = opts, target_mean = 0.71, target_x = c(0.1, 0.9), target_y = c(0.05, 0.95))

mlog(c(0.1, 0.9), 0.90, mean = 0.71, n = 3, bounds = c(0, 1), method = "nl")

res
res$obj
ml = mlog(res$sol, "coeff", bounds = c(0, 1))
plot(ml)
moment.mlog(ml)
qmlog(0.05, ml)
qmlog(0.95, ml)

ml = mlog(c(0, 1, 0), "coeff", bounds = c(0, 1))
ml$a = res$sol
plot(qmlog(seq(0, 1, by = 0.01), ml))
