# Code for metalogs specified via symmetric-percentile triplets, i.e. for
# y = c(p, 1/2, 1 - p)
SPT_a = function(q1, q2, q3, p, bounds)
{
    bl = bounds[1];
    bu = bounds[2];

    if (is.finite(bl)) {
        if (is.finite(bu)) {
            # c(bl, bu)
            # Proposition 3, Keelin 2016, p.270:
            q1 = (q1 - bl) / (bu - q1);
            q2 = (q2 - bl) / (bu - q2);
            q3 = (q3 - bl) / (bu - q3);
            a = c(
                log(q2),
                1 / (2 * log((1-p)/p)) * log(q3/q1),
                log(q3 * q1 / q2^2) / ((1 - 2*p) * log((1-p)/p))
            );
        } else {
            # c(bl, Inf)
            # Proposition 3, Keelin 2016, p.270:
            q1 = q1 - bl;
            q2 = q2 - bl;
            q3 = q3 - bl;
            a = c(
                log(q2),
                1 / (2 * log((1-p)/p)) * log(q3/q1),
                log(q3 * q1 / q2^2) / ((1 - 2*p) * log((1-p)/p))
            );
        }
    } else {
        if (is.finite(bu)) {
            # c(-Inf, bu)
            q1 = bu - q1;
            q2 = bu - q2;
            q3 = bu - q3;
            a = c(
                -log(q2),
                -1 / (2 * log((1-p)/p)) * log(q3/q1),
                -log(q3 * q1 / q2^2) / ((1 - 2*p) * log((1-p)/p))
            );
        } else {
            # c(-Inf, Inf)
            # eq 16, Keelin 2016:
            a = c(
                q2,
                1 / (2 * log((1-p)/p)) * (q3 - q1),
                (q3 + q1 - 2 * q2) / ((1 - 2*p) * log((1-p)/p))
            )
        }
    }

    return (a)
}
