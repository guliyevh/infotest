#' Information Matrix (IM) Test for \code{lm()} Models
#'
#' @title Information Matrix Test for Regression Models
#'
#' @description
#' \code{infotest()} applies the Information Matrix (IM) test to a fitted
#' \code{\link[stats]{lm}} model and reports the Cameron & Trivedi (1990) decomposition
#' into heteroskedasticity, skewness, and kurtosis components. Optionally,
#' White's (1980) heteroskedasticity test is also reported.
#'
#' @param model A fitted \code{\link[stats]{lm}} object.
#' @param white Logical; if \code{TRUE}, also compute and report White's test
#'   (default \code{TRUE}).
#'
#' @return
#' Invisibly returns a \code{list} with:
#' \itemize{
#'   \item \code{decomposition}: a list with components
#'     \code{heteroskedasticity}, \code{skewness}, \code{kurtosis}, and
#'     \code{total}, each providing \code{chi2}, \code{df}, and \code{p}.
#'   \item \code{white}: (only when \code{white = TRUE}) a list with
#'     \code{chi2}, \code{df}, and \code{p}.
#' }
#'
#' @details
#' The function checks that \code{model} is an \code{lm} fit and refuses
#' weighted models. The intercept is removed from the design matrix, regressors
#' are centered, and quadratic terms are constructed for the auxiliary
#' regressions used by the IM and White tests. Test statistics are computed as
#' chi-square values with associated degrees of freedom. A formatted table is
#' printed to the console; the full results are returned invisibly.
#'
#' @section Warnings:
#' \itemize{
#'   \item Weighted \code{lm} fits are not supported.
#'   \item If there are no covariates (intercept-only model), the test is
#'     skipped with a message.
#' }
#'
#' @references
#' White, H. (1980). A heteroskedasticity-consistent covariance matrix
#' estimator and a direct test for heteroskedasticity. \emph{Econometrica},
#' 48(4), 817–838. \cr
#' Cameron, A. C., & Trivedi, P. K. (1990). The information matrix test and its
#' applied alternatives. \emph{Econometric Theory}, 6(2), 179–195.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{pchisq}}
#'
#' @keywords htest diagnostics
#' @concept diagnostics
#' @author Hasraddin Guliyev \email{hasradding@unec.edu.az}
#' @rdname infotest
#'
#' @examples
#' m <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)
#' infotest(m)
#' infotest(m, white = TRUE)
#' res <- invisible(infotest(m))
#' str(res$decomposition)
#'
#' @importFrom stats lm residuals model.response model.frame model.matrix pchisq
#' @export
infotest <- function(model, white = TRUE) {
  # Check if the model is from lm()
  if (!inherits(model, "lm")) {
    stop("imtest is only possible after linear regression (lm)")
  }

  # Check for weights (not supported)
  if (!is.null(model$weights)) {
    stop("imtest does not support weights")
  }

  # Extract model components
  y <- model.response(model.frame(model))
  X <- model.matrix(model)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE] # Remove intercept
  nrhs <- ncol(X)

  if (nrhs == 0) {
    message("(imtest not allowed without covariates)")
    return(invisible(NULL))
  }

  res <- residuals(model)
  n <- length(res)
  s2 <- sum(res^2) / n

  # Create squared terms (x_i * x_j)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  rhs2 <- NULL
  k <- 0
  for (i in 1:nrhs) {
    for (j in 1:i) {
      k <- k + 1
      rhs2 <- cbind(rhs2, X_centered[, i] * X_centered[, j])
    }
  }

  results <- list()

  # White's original test
  if (white) {
    y_white <- res^2
    white_model <- lm(y_white ~ X + rhs2)
    r2 <- summary(white_model)$r.squared
    chi2 <- n * r2
    df <- ncol(X) + ncol(rhs2)
    p <- pchisq(chi2, df, lower.tail = FALSE)

    results$white <- list(
      chi2 = chi2,
      df = df,
      p = p
    )
  }

  # Cameron & Trivedi's decomposition

  # Heteroskedasticity component
  y_h <- res^2 - s2
  y_h_sq <- y_h^2
  uss_h <- sum(y_h_sq)
  model_h <- lm(y_h ~ X + rhs2)
  rss_h <- sum(residuals(model_h)^2)
  chi2_h <- n * (1 - rss_h / uss_h)
  df_h <- model_h$rank - 1
  p_h <- pchisq(chi2_h, df_h, lower.tail = FALSE)

  # Skewness component
  y_s <- res^3 - 3 * s2 * res
  y_s_sq <- y_s^2
  uss_s <- sum(y_s_sq)
  model_s <- lm(y_s ~ X)
  rss_s <- sum(residuals(model_s)^2)
  chi2_s <- n * (1 - rss_s / uss_s)
  df_s <- model_s$rank - 1
  p_s <- pchisq(chi2_s, df_s, lower.tail = FALSE)

  # Kurtosis component
  y_k <- res^4 - 6 * s2 * res^2 + 3 * s2^2
  y_k_sq <- y_k^2
  uss_k <- sum(y_k_sq)
  model_k <- lm(y_k ~ 1)
  rss_k <- sum(residuals(model_k)^2)
  chi2_k <- n * (1 - rss_k / uss_k)
  df_k <- 1
  p_k <- pchisq(chi2_k, df_k, lower.tail = FALSE)

  # Total
  chi2_t <- chi2_s + chi2_k + chi2_h
  df_t <- df_s + df_k + df_h
  p_t <- pchisq(chi2_t, df_t, lower.tail = FALSE)

  # Store results
  results$decomposition <- list(
    heteroskedasticity = list(chi2 = chi2_h, df = df_h, p = p_h),
    skewness = list(chi2 = chi2_s, df = df_s, p = p_s),
    kurtosis = list(chi2 = chi2_k, df = df_k, p = p_k),
    total = list(chi2 = chi2_t, df = df_t, p = p_t)
  )

  # Print results
  if (white) {
    cat("\nWhite's (1980) heteroskedasticity test\n")
    cat("H0: Homoskedasticity\n")
    cat("Ha: Unrestricted heteroskedasticity\n\n")
    cat(sprintf("chi2(%d) = %.4f\n", results$white$df, results$white$chi2))
    cat(sprintf("Prob > chi2 = %.4f\n\n", results$white$p))
  }

  cat("Cameron & Trivedi's (1990) decomposition of IM-test\n\n")

  # Create table
  tab <- data.frame(
    Source = c("Heteroskedasticity", "Skewness", "Kurtosis", "Joint Test"),
    chi2 = c(chi2_h, chi2_s, chi2_k, chi2_t),
    df = c(df_h, df_s, df_k, df_t),
    p = c(p_h, p_s, p_k, p_t)
  )

  # Print formatted table
  cat("----------------------------------------\n")
  cat(sprintf("%-20s %10s %6s %10s\n", "Source", "chi2", "df", "p"))
  cat("----------------------------------------\n")
  for (i in 1:nrow(tab)) {
    cat(sprintf("%-20s %10.2f %6d %10.4f\n",
                tab$Source[i], tab$chi2[i], tab$df[i], tab$p[i]))
  }
  cat("----------------------------------------\n")

  return(invisible(results))
}
