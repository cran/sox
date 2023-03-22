#' plots for \code{sox_cv}
#'
#' @description
#' Plot the cross-validation curves produced by \code{\link{sox_cv}}.
#' 
#' @param sox_cv_obj The \code{\link{sox_cv}} object.
#' 
#' @return
#'The plot is the \code{cvm} (red dot) for each lambda with its standard error (vertical bar). The two vertical dashed lines corresponds to the \code{lambda.min} and \code{lambda.1se}
#' 
#' @examples 
#' grp <- matrix(c(0, 0, 0, 0, 0,
#'                 0, 0, 0, 0, 0,
#'                 1, 1, 0, 0, 0,
#'                 0, 0, 0, 0, 0,
#'                 0, 1, 0, 1, 0),
#'               ncol = 5, byrow = TRUE)
#' grp.var <- matrix(c(1, 0, 0, 0, 0, #A1
#'                     1, 0, 0, 0, 0, #A2
#'                     0, 0, 0, 1, 0, #C1
#'                     0, 0, 0, 1, 0, #C2
#'                     0, 1, 0, 0, 0, #B
#'                     0, 0, 1, 0, 0, #A1B
#'                     0, 0, 1, 0, 0, #A2B
#'                     0, 0, 0, 0, 1, #C1B
#'                     0, 0, 0, 0, 1  #C2B
#'                    ), ncol = 5, byrow = TRUE)
#' eta_g <- rep(1, 5)
#' x <- as.matrix(sim[, c("A1","A2","C1","C2","B",
#'                        "A1B","A2B","C1B","C2B")])
#' lam.seq <- 10^seq(0, -2, by = -0.2)
#' 
#' cv <- sox_cv(x = x,
#'                 ID = sim$Id,
#'                 time = sim$Start,
#'                 time2 = sim$Stop,
#'                 event = sim$Event,
#'                 lambda = lam.seq,
#'                 group = grp,
#'                 group_variable = grp.var,
#'                 penalty_weights = eta_g,
#'                 nfolds = 5,
#'                 tol = 1e-4,
#'                 maxit = 1e3,
#'                 verbose = FALSE)
#' plot_sox_cv(cv)               
#' @seealso \code{\link{sox}}, \code{\link{sox_cv}}.

plot_sox_cv <- function(sox_cv_obj) {
  cvup <- sox_cv_obj$cvup
  cvlo <- sox_cv_obj$cvlo
  cvm <- sox_cv_obj$cvm
  
  lam.seq <- sox_cv_obj$lambdas
  lambda.min <- sox_cv_obj$lambda.min
  lambda.1se <- sox_cv_obj$lambda.1se
  
  ylims <- range(c(cvup, cvlo))
  plot(x = log(lam.seq), y = cvm, ylim = ylims, col = 0,
       ylab = "CV Error",
       xlab = expression(paste("Log(", lambda, ")")))
  for (i in 1:length(cvm)) {
    segments(x0 = log(lam.seq)[i], y0 = cvup[i],
             x1 = log(lam.seq)[i], y1 = cvlo[i],
             col = 'grey60')
  }
  points(x = log(lam.seq), y = cvup, pch = 95, cex = 1, col = 'grey60')
  points(x = log(lam.seq), y = cvlo, pch = 95, cex = 1, col = 'grey60')
  points(x = log(lam.seq), y = cvm, pch = 20, cex = 1, col = 'red')
  abline(v = log(lambda.min), lty = 3)
  abline(v = log(lambda.1se), lty = 3)
  axis(side = 3, at = log(lam.seq), labels = sox_cv_obj$nzero, tick = FALSE)
}
