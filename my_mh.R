my_mantelhaen.test <-
function (x, y = NULL, z = NULL, alternative = c("two.sided", 
    "less", "greater"), correct = TRUE, exact = FALSE, conf.level = 0.95) 
{
    DNAME <- deparse(substitute(x))
    if (is.array(x)) {
        if (length(dim(x)) == 3L) {
            if (anyNA(x)) 
                stop("NAs are not allowed")
            if (any(dim(x) < 2L)) 
                stop("each dimension in table must be >= 2")
        }
        else stop("'x' must be a 3-dimensional array")
    }
    else {
        if (is.null(y)) 
            stop("if 'x' is not an array, 'y' must be given")
        if (is.null(z)) 
            stop("if 'x' is not an array, 'z' must be given")
        if (any(diff(c(length(x), length(y), length(z))) != 0L)) 
            stop("'x', 'y', and 'z' must have the same length")
        DNAME <- paste(DNAME, "and", deparse(substitute(y)), 
            "and", deparse(substitute(z)))
        OK <- complete.cases(x, y, z)
        x <- factor(x[OK])
        y <- factor(y[OK])
        if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
            stop("'x' and 'y' must have at least 2 levels")
        else x <- table(x, y, z[OK])
    }
    if (any(apply(x, 3L, sum) < 2)) 
        stop("sample size in each stratum must be > 1")
    I <- dim(x)[1L]
    J <- dim(x)[2L]
    K <- dim(x)[3L]
    if ((I == 2) && (J == 2)) {
        alternative <- match.arg(alternative)
        if (!missing(conf.level) && (length(conf.level) != 1 || 
            !is.finite(conf.level) || conf.level < 0 || conf.level > 
            1)) 
            stop("'conf.level' must be a single number between 0 and 1")
        NVAL <- c(`common odds ratio` = 1)
        if (!exact) {
            s.x <- apply(x, c(1L, 3L), sum)
            s.y <- apply(x, c(2L, 3L), sum)
            n <- as.double(apply(x, 3L, sum))
            DELTA <- sum(x[1, 1, ] - s.x[1, ] * s.y[1, ]/n)
            YATES <- if (correct && (abs(DELTA) >= 0.5)) 
                0.5
            else 0
            STATISTIC <- ((abs(DELTA) - YATES)^2/sum(apply(rbind(s.x, 
                s.y), 2L, prod)/(n^2 * (n - 1))))
            PARAMETER <- 1
            if (alternative == "two.sided") 
                PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
            else {
                z <- sign(DELTA) * sqrt(STATISTIC)
                PVAL <- pnorm(z, lower.tail = (alternative == 
                  "less"))
            }
            names(STATISTIC) <- "Mantel-Haenszel X-squared"
            names(PARAMETER) <- "df"
            METHOD <- paste("Mantel-Haenszel chi-squared test", 
                if (YATES) 
                  "with"
                else "without", "continuity correction")
            s.diag <- sum(x[1L, 1L, ] * x[2L, 2L, ]/n)
            s.offd <- sum(x[1L, 2L, ] * x[2L, 1L, ]/n)
            ESTIMATE <- s.diag/s.offd
            sd <- sqrt(sum((x[1L, 1L, ] + x[2L, 2L, ]) * x[1L, 
                1L, ] * x[2L, 2L, ]/n^2)/(2 * s.diag^2) + sum(((x[1L, 
                1L, ] + x[2L, 2L, ]) * x[1L, 2L, ] * x[2L, 1L, 
                ] + (x[1L, 2L, ] + x[2L, 1L, ]) * x[1L, 1L, ] * 
                x[2L, 2L, ])/n^2)/(2 * s.diag * s.offd) + sum((x[1L, 
                2L, ] + x[2L, 1L, ]) * x[1L, 2L, ] * x[2L, 1L, 
                ]/n^2)/(2 * s.offd^2))
            CINT <- switch(alternative, less = c(0, ESTIMATE * 
                exp(qnorm(conf.level) * sd)), greater = c(ESTIMATE * 
                exp(qnorm(conf.level, lower.tail = FALSE) * sd), 
                Inf), two.sided = {
                ESTIMATE * exp(c(1, -1) * qnorm((1 - conf.level)/2) * 
                  sd)
            })
            RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
                p.value = PVAL)
        }
        else {
            METHOD <- paste("Exact conditional test of independence", 
                "in 2 x 2 x k tables")
            mn <- apply(x, c(2L, 3L), sum)
            m <- mn[1L, ]
            n <- mn[2L, ]
            t <- apply(x, c(1L, 3L), sum)[1L, ]
            s <- sum(x[1L, 1L, ])
            lo <- sum(pmax(0, t - n))
            hi <- sum(pmin(m, t))
            support <- lo:hi
            dc <- .Call(C_d2x2xk, K, m, n, t, hi - lo + 1L)
            logdc <- log(dc)
            dn2x2xk <- function(ncp) {
                if (ncp == 1) 
                  return(dc)
                d <- logdc + log(ncp) * support
                d <- exp(d - max(d))
                d/sum(d)
            }
            mn2x2xk <- function(ncp) {
                if (ncp == 0) 
                  return(lo)
                if (ncp == Inf) 
                  return(hi)
                sum(support * dn2x2xk(ncp))
            }
            pn2x2xk <- function(q, ncp = 1, upper.tail = FALSE) {
                if (ncp == 0) {
                  if (upper.tail) 
                    return(as.numeric(q <= lo))
                  else return(as.numeric(q >= lo))
                }
                if (ncp == Inf) {
                  if (upper.tail) 
                    return(as.numeric(q <= hi))
                  else return(as.numeric(q >= hi))
                }
                d <- dn2x2xk(ncp)
                if (upper.tail) 
                  sum(d[support >= q])
                else sum(d[support <= q])
            }
            PVAL <- switch(alternative, less = pn2x2xk(s, 1), 
                greater = pn2x2xk(s, 1, upper.tail = TRUE), two.sided = {
                  relErr <- 1 + 10^(-7)
                  d <- dc
                  sum(d[d <= d[s - lo + 1] * relErr])
                })
            mle <- function(x) {
                if (x == lo) 
                  return(0)
                if (x == hi) 
                  return(Inf)
                mu <- mn2x2xk(1)
                if (mu > x) 
                  uniroot(function(t) mn2x2xk(t) - x, c(0, 1))$root
                else if (mu < x) 
                  1/uniroot(function(t) mn2x2xk(1/t) - x, c(.Machine$double.eps, 
                    1))$root
                else 1
            }
            ESTIMATE <- mle(s)
            ncp.U <- function(x, alpha) {
                if (x == hi) 
                  return(Inf)
                p <- pn2x2xk(x, 1)
                if (p < alpha) 
                  uniroot(function(t) pn2x2xk(x, t) - alpha, 
                    c(0, 1))$root
                else if (p > alpha) 
                  1/uniroot(function(t) pn2x2xk(x, 1/t) - alpha, 
                    c(.Machine$double.eps, 1))$root
                else 1
            }
            ncp.L <- function(x, alpha) {
                if (x == lo) 
                  return(0)
                p <- pn2x2xk(x, 1, upper.tail = TRUE)
                if (p > alpha) 
                  uniroot(function(t) pn2x2xk(x, t, upper.tail = TRUE) - 
                    alpha, c(0, 1))$root
                else if (p < alpha) 
                  1/uniroot(function(t) pn2x2xk(x, 1/t, upper.tail = TRUE) - 
                    alpha, c(.Machine$double.eps, 1))$root
                else 1
            }
            CINT <- switch(alternative, less = c(0, ncp.U(s, 
                1 - conf.level)), greater = c(ncp.L(s, 1 - conf.level), 
                Inf), two.sided = {
                alpha <- (1 - conf.level)/2
                c(ncp.L(s, alpha), ncp.U(s, alpha))
            })
            STATISTIC <- c(S = s)
            RVAL <- list(statistic = STATISTIC, p.value = PVAL)
        }
        names(ESTIMATE) <- names(NVAL)
        attr(CINT, "conf.level") <- conf.level
        RVAL <- c(RVAL, list(conf.int = CINT, estimate = ESTIMATE, 
            null.value = NVAL, alternative = alternative))
    }
    else {
        df <- (I - 1) * (J - 1)
        n <- m <- double(length = df)
        V <- matrix(0, nrow = df, ncol = df)
        for (k in 1:K) {
            f <- x[, , k]
            ntot <- sum(f)
            rowsums <- apply(f, 1L, sum)[-I]
            colsums <- apply(f, 2L, sum)[-J]
            n <- n + c(f[-I, -J])
            ntot <- as.numeric(ntot) # additional info
            m <- m + c(outer(rowsums, colsums, "*"))/ntot
            V <- V + (kronecker(diag(ntot * colsums, nrow = J - 
                1) - outer(colsums, colsums), diag(ntot * rowsums, 
                nrow = I - 1) - outer(rowsums, rowsums))/(ntot^2 * 
                (ntot - 1)))
        }
        n <- n - m
        STATISTIC <- c(crossprod(n, qr.solve(V, n)))
        PARAMETER <- df
        PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
        names(STATISTIC) <- "Cochran-Mantel-Haenszel M^2"
        names(PARAMETER) <- "df"
        METHOD <- "Cochran-Mantel-Haenszel test"
        RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
            p.value = PVAL)
    }
    RVAL <- c(RVAL, list(method = METHOD, data.name = DNAME))
    class(RVAL) <- "htest"
    return(RVAL)
}
