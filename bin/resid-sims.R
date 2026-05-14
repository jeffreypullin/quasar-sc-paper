# nolint start

N <- 1000
n <- 20
w <- rnorm(N, 10)
X <- cbind(rep(1, N), rnorm(N), rnorm(N))
sample_id <- as.factor(rep(1:n, each = N / n))
W <- diag(w)
Z <- model.matrix(~ 0 + sample_id)

PW_full <- diag(N) - X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

XtWX_inv <- solve(t(X) %*% W %*% X)
W_summed <- t(Z) %*% W %*% Z
X_mean <- (t(Z) %*% X) / (N / n)
X_summed <-  (t(Z) %*% X)
X_w_mean <-  (t(Z) %*% W %*% X) / (diag(t(Z) %*% W %*% Z))

XtWX_inv
solve(t(X_mean) %*% W_summed %*% X_mean)
solve(t(X_summed) %*% W_summed %*% X_summed)

test1 <- t(Z) %*% PW_full
test2 <- (diag(n) - X_mean %*% XtWX_inv %*% t(X_mean) %*% W_summed) %*% t(Z)
test2 <- unname(test2)
test3 <- (diag(n) - X_mean %*% solve(t(X_mean) %*% W_summed %*% X_mean) %*% t(X_mean) %*% W_summed) %*% t(Z)
test3 <- unname(test3)
test4 <- (diag(n) - X_mean %*% solve(t(X_summed) %*% W_summed %*% X_summed) %*% t(X_mean) %*% W_summed) %*% t(Z)
test4 <- unname(test4)
test5 <- (diag(n) - X_w_mean %*% solve(t(X_w_mean) %*% W_summed %*% X_w_mean) %*% t(X_w_mean) %*% W_summed) %*% t(Z)
test5 <- unname(test5)

W %*% Z


test1[1:5, 1:5]
test2[1:5, 1:5]
test3[1:5, 1:5]
test4[1:5, 1:5]
test5[1:5, 1:5]

(t(Z) %*% t(test1))[1:5, 1:5]
test2[1:5, 1:5]

as.vector((t(Z) %*% test1[1, ]) / 50)
unique(test2[1, ])

((t(Z) %*% PW_full %*% Z) / 50)[1:5, 1:5]
test2[1:5, 1:5]

# nolint end
