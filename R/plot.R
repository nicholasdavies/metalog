#' @export
plot.mlog = function(ml)
{
    saved_mfrow = par("mfrow");
    saved_mar = par("mar");
    on.exit(par(mfrow = saved_mfrow, mar = saved_mar));

    par(mfrow = c(2, 1));
    par(mar = c(2.1, 4.1, 2.1, 2.1));

    M = ml$cache$M;
    m = ml$cache$m;
    y = ml$cache$y;

    plot(M, y, type = "l", main = "CDF", xlab = NA, ylab = "Probability", ylim = c(0, 1));
    polygon(
        c(M[1], M, M[length(M)]),
        c(0, y, 0),
        col = rgb(0.75, 0.75, 0.75), border = NA
    );

    plot(M, m, type = "l", main = "PDF", xlab = NA, ylab = "Density", ylim = c(0, max(m)));
    polygon(
        c(M[1], M, M[length(M)]),
        c(0, m, 0),
        col = rgb(0.75, 0.75, 0.75), border = NA
    );
}
