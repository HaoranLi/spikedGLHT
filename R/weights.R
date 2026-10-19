# Corollary 3.1 and Proposition 4.1.
.weight_geometry <- function(ell, sigma_sq, n, p) {
    K <- length(ell)
    z2 <- pmax(0, (1 - (p / n) / (ell / sigma_sq - 1)^2) /
                   (1 + (p / n) / (ell / sigma_sq - 1)))
    z2[z2 < 1e-12] <- 0
    a <- (z2 * ell^2 + (1 - z2) * sigma_sq^2) / p
    d <- (z2 * ell + (1 - z2) * sigma_sq)^2 / p
    G <- diag(c((sum(ell^2) + (p - K) * sigma_sq^2) / p, d), K + 1L)
    if (K > 0L) {
        G[1L, -1L] <- a
        G[-1L, 1L] <- a
    }
    list(G = G, z2 = z2)
}

.minimax <- function(geometry, rho) {
    G <- geometry$G
    z2 <- geometry$z2
    K <- length(z2)
    if (K == 0L) return(list(weights = 1, prior = 1))
    active <- which(z2 > 0)
    inactive <- which(z2 == 0)
    # At the detection boundary, equal b columns are combined exactly.
    # Their mass can be assigned to an undetectable spike with pi_0=0.
    Z <- diag(c(1, z2[active]), length(active) + 1L)
    Z[1L, ] <- 1
    rows <- c(1L, active + 1L)
    H <- crossprod(Z, solve(G)[rows, rows, drop = FALSE] %*% Z)
    size <- nrow(H)
    constraints <- cbind(rep(1, size), diag(size))
    bounds <- c(1, rep(0, size))
    if (!length(inactive)) {
        constraints <- cbind(constraints, c(-1, rep(0, size - 1L)))
        bounds <- c(bounds, -rho)
    }
    # Rescaling preserves the minimizer and improves QP conditioning.
    H <- H / max(abs(H))
    solution <- quadprog::solve.QP(H, rep(0, size), constraints, bounds, meq = 1)$solution
    prior <- numeric(K + 1L)
    prior[active + 1L] <- solution[-1L]
    prior[if (length(inactive)) inactive[1L] + 1L else 1L] <- solution[1L]
    b <- c(sum(prior), z2 * prior[-1L])
    weights <- drop(solve(G, b))
    weights <- weights / sqrt(sum(weights^2))
    list(weights = weights, prior = prior)
}
