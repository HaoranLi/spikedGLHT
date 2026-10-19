.check_scalar <- function(x, name, lower, upper, integer = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > upper || (integer && x != floor(x)))
        stop(name, " must be a finite ", if (integer) "integer " else "number ",
             "between ", lower, " and ", upper, ".", call. = FALSE)
}

.check_matrix <- function(x, name) {
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
        any(dim(x) == 0L))
        stop(name, " must be a nonempty finite numeric matrix.", call. = FALSE)
}

.inverse_sqrt <- function(x) {
    e <- eigen(x, symmetric = TRUE)
    if (min(e$values) <= max(e$values) * .Machine$double.eps * nrow(x))
        stop("The contrast covariance is numerically singular.", call. = FALSE)
    tcrossprod(sweep(e$vectors, 2L, sqrt(e$values), `/`), e$vectors)
}

# Algorithm 1. The optional N/n factor reproduces the research scripts.
.adjusted <- function(values, K, N, n, p, adjustment) {
    if (K == 0L) return(list(spikes = numeric(), sigma_sq = mean(values)))
    bulk <- mean(values[seq.int(K + 1L, p)])
    if (bulk <= 0) stop("Residual bulk variance must be positive.", call. = FALSE)
    gamma <- p / n
    invert <- function(sigma_sq, final = FALSE) {
        z <- values[seq_len(K)] / sigma_sq
        if (adjustment == "legacy") z <- z * N / n
        b <- gamma - 1 - z
        root <- (-b + sqrt(pmax(b^2 - 4 * z, 0))) / 2
        if (adjustment == "legacy") {
            root <- pmax(root, 1 + sqrt(gamma) + 1e-4)
        } else if (final) {
            root[z < (1 + sqrt(gamma))^2] <- 1 + sqrt(gamma)
        }
        root * sigma_sq
    }
    initial <- invert(bulk)
    b0 <- K + sum(1 / (initial / bulk - 1))
    sigma_sq <- bulk * (1 + gamma * b0 / (p - K))
    if (!is.finite(sigma_sq) || sigma_sq <= 0)
        stop("Noise estimation failed; reduce K or inspect the residual spectrum.", call. = FALSE)
    list(spikes = invert(sigma_sq, TRUE), sigma_sq = sigma_sq)
}

.select_K <- function(values, N, n, p, method, conf_level, delta, K_max, adjustment) {
    if (K_max == 0L) return(0L)
    noise <- function(k) .adjusted(values, k, N, n, p, adjustment)$sigma_sq
    if (method == "PLY") {
        penalty <- (p / n + 2 * sqrt(p / n)) * (1 + n / p^(1 + delta)) / p
        scores <- vapply(0:K_max, noise, numeric(1)) +
            (0:K_max) * noise(K_max) * penalty
        return(which.min(scores) - 1L)
    }
    if (!requireNamespace("RMTstat", quietly = TRUE))
        stop("Install RMTstat for KN spike selection, or supply K or use PLY.", call. = FALSE)
    s <- RMTstat::qtw(conf_level, beta = 1)
    # Preserve the research implementation's sequential test, while allowing K=0.
    for (k in seq_len(K_max)) {
        d <- p - k
        mu <- (sqrt(n - 0.5) + sqrt(d - 0.5))^2 / n
        eta <- (sqrt(n - 0.5) + sqrt(d - 0.5)) *
            (1 / sqrt(n - 0.5) + 1 / sqrt(d - 0.5))^(1 / 3) / n
        if (values[k] <= noise(k) * (mu + eta * s)) return(k - 1L)
    }
    K_max
}
