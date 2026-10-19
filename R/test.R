.raw_statistics <- function(M) {
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev <= -1))
        stop("LR and BNP statistics are outside their domain (an eigenvalue is <= -1).",
             call. = FALSE)
    c(LR = sum(log1p(ev)), LH = sum(ev), BNP = sum(ev / (1 + ev)))
}

#' Test a general linear hypothesis with spiked noise covariance
#' @export
spiked_glht <- function(Y, X, C, D = NULL, rho = 0.5, K = NULL,
                        calibration = c("normal", "bootstrap"), alpha = 0.05,
                        B = 1999L, method = c("KN", "PLY"),
                        conf_level = 0.95, delta = 0.05, K_max = NULL) {
    calibration <- match.arg(calibration)
    method <- match.arg(method)
    .check_matrix(Y, "Y"); .check_matrix(X, "X"); .check_matrix(C, "C")
    p <- nrow(Y); N <- ncol(Y); m <- nrow(X); n <- N - m; q <- ncol(C)
    if (ncol(X) != N || nrow(C) != m || m >= N || q > m)
        stop("Require Y: p x N, X: m x N, C: m x q, q <= m < N.", call. = FALSE)
    .check_scalar(rho, "rho", .Machine$double.eps, 1)
    .check_scalar(alpha, "alpha", .Machine$double.eps, 1 - .Machine$double.eps)
    .check_scalar(conf_level, "conf_level", .Machine$double.eps, 1 - .Machine$double.eps)
    .check_scalar(delta, "delta", .Machine$double.eps, Inf)
    if (calibration == "bootstrap") .check_scalar(B, "B", 1, .Machine$integer.max, TRUE)
    qrX <- qr(t(X))
    if (qrX$rank != m || qr(C)$rank != q)
        stop("X must have full row rank and C full column rank.", call. = FALSE)
    if (is.null(D)) D <- matrix(0, p, q)
    .check_matrix(D, "D")
    if (!identical(dim(D), c(p, q))) stop("D must have dimensions p x q.", call. = FALSE)
    R <- .inverse_sqrt(crossprod(C, solve(tcrossprod(X), C)))
    # QR residuals avoid forming an N x N projection. SVD also avoids p x p covariance.
    residual <- qr.resid(qrX, t(Y))
    spectrum <- svd(residual / sqrt(n), nu = 0L, nv = min(p, n))
    values <- c(spectrum$d^2, rep(0, max(0, p - length(spectrum$d))))
    limit <- min(p - 1L, n - 1L)
    if (is.null(K_max)) K_max <- min(limit, max(1L, floor(p * 0.1)))
    .check_scalar(K_max, "K_max", 0, limit, TRUE)
    if (is.null(K)) K <- .select_K(values, n, p, method, conf_level,
                                  delta, K_max)
    .check_scalar(K, "K", 0, limit, TRUE)
    if (calibration == "bootstrap" && p < 2 * K)
        stop("Bootstrap representation requires p >= 2*K.", call. = FALSE)
    if (K > 0L && values[K] <= max(values) * .Machine$double.eps)
        stop("K exceeds the numerical rank of the residuals.", call. = FALSE)
    fit <- .adjusted(values, K, n, p)
    if (fit$sigma_sq <= 0) stop("Residual variance must be positive.", call. = FALSE)
    geometry <- .weight_geometry(fit$spikes, fit$sigma_sq, n, p)
    opt <- .minimax(geometry, rho)
    lambda <- opt$weights
    theta1 <- (sum(lambda[-1L] * (geometry$z2 * fit$spikes +
                  (1 - geometry$z2) * fit$sigma_sq)) +
                  lambda[1L] * (sum(fit$spikes) + (p - K) * fit$sigma_sq)) / p
    theta2 <- drop(crossprod(lambda, geometry$G %*% lambda))
    Q <- t(X) %*% solve(tcrossprod(X), C) %*% R
    A <- Y %*% Q - D %*% R
    M <- lambda[1L] * crossprod(A)
    if (K > 0L) {
        projected <- crossprod(spectrum$v[, seq_len(K), drop = FALSE], A)
        M <- M + crossprod(projected, sweep(projected, 1L, lambda[-1L], `*`))
    }
    raw <- .raw_statistics(M / p)
    if (theta2 <= 0 || 1 + theta1 <= 0)
        stop("Invalid estimated centering or variance.", call. = FALSE)
    center <- c(LR = q * log1p(theta1), LH = q * theta1,
                BNP = q * theta1 / (1 + theta1))
    scale <- sqrt(p / (2 * q * theta2)) *
        c(LR = 1 + theta1, LH = 1, BNP = (1 + theta1)^2)
    statistic <- (raw - center) * scale
    normal_p <- stats::pnorm(statistic, lower.tail = FALSE)
    boot <- NULL
    p_value <- normal_p
    critical <- center + stats::qnorm(1 - alpha) / scale
    if (calibration == "bootstrap") {
        boot <- .bootstrap(B, q, p, fit, geometry$z2, lambda)
        p_value <- (1 + colSums(sweep(boot, 2L, raw, `>=`))) / (B + 1)
        critical <- apply(boot, 2L, stats::quantile, probs = 1 - alpha, names = FALSE)
    }
    names(lambda) <- names(opt$prior) <- paste0("component", 0:K)
    structure(list(statistic = statistic, p.value = p_value,
                   reject = raw > critical, raw.statistic = raw,
                   critical.value = critical, normal.p.value = normal_p,
                   weights = lambda, least.favorable.prior = opt$prior,
                   K = K, spikes = fit$spikes, sigma_sq = fit$sigma_sq,
                   theta1 = theta1, theta2 = theta2, alpha = alpha,
                   calibration = calibration,
                   rho = rho, dimensions = c(p = p, N = N, m = m, n = n, q = q),
                   bootstrap.statistics = boot, call = match.call()), class = "spiked_glht")
}

.bootstrap <- function(B, q, p, fit, z2, lambda) {
    K <- length(z2)
    ell <- fit$spikes; sigma_sq <- fit$sigma_sq
    a1 <- (lambda[1L] + lambda[-1L] * z2) * ell / p
    a2 <- (lambda[1L] + lambda[-1L] * (1 - z2)) * sigma_sq / p
    a3 <- sqrt(sigma_sq * ell) * lambda[-1L] * sqrt(z2 * (1 - z2)) / p
    out <- matrix(0, B, 3L, dimnames = list(NULL, c("LR", "LH", "BNP")))
    for (i in seq_len(B)) {
        V1 <- matrix(stats::rnorm(K * q), K, q)
        V2 <- matrix(stats::rnorm(K * q), K, q)
        V3 <- matrix(stats::rnorm((p - 2 * K) * q), p - 2 * K, q)
        M <- sigma_sq * lambda[1L] / p * crossprod(V3)
        if (K > 0L) {
            cross <- crossprod(V1, sweep(V2, 1L, a3, `*`))
            M <- M + crossprod(V1, sweep(V1, 1L, a1, `*`)) +
                crossprod(V2, sweep(V2, 1L, a2, `*`)) + cross + t(cross)
        }
        out[i, ] <- .raw_statistics(M)
    }
    out
}

#' @export
print.spiked_glht <- function(x, ...) {
    cat("General linear hypothesis test with spiked noise covariance\n")
    cat("Calibration:", x$calibration, " | K:", x$K, " | rho:", x$rho, "\n")
    print(data.frame(statistic = x$statistic, p.value = x$p.value, reject = x$reject), ...)
    invisible(x)
}
