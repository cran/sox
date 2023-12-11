#' Plots for \code{sox_cv}
#'
#' @description
#' Plot the solution path or cross-validation curves produced by \code{\link{sox_cv}()}.
#' 
#' @param x The \code{\link{sox_cv}} object.
#' @param type Character string, "\code{solution-path}" to generate a solution path with marks at \code{lambda.min} and \code{lambda.1se}; "\code{cv-curve}" to generate a cross-validation curve.
#' @param ... Other graphical parameters to plot
#' 
#' @return
#' The "\code{solution-path}" plot produces a coefficient profile plot of the coefficient paths for a fitted \code{\link{sox}} model. The "\code{cv-curve}" plot is the \code{cvm} (red dot) for each lambda with its standard error (vertical bar). The two vertical dashed lines corresponds to the \code{lambda.min} and \code{lambda.1se}
#' 
#' @examples 
#' x <- as.matrix(sim[, c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B")])
#' lam.seq <- exp(seq(log(1e0), log(1e-3), length.out = 20))
#' 
#' overlapping.groups <- list(c(1, 2, 6, 7),
#'                            c(5, 6, 7, 8, 9),
#'                            c(6, 7),
#'                            c(3, 4, 8, 9),
#'                            c(8, 9))
#'                            
#' pars.overlapping <- overlap_structure(overlapping.groups)
#'               
#' cv.overlapping <- sox_cv(
#'   x = x,
#'   ID = sim$Id,
#'   time = sim$Start,
#'   time2 = sim$Stop,
#'   event = sim$Event,
#'   penalty = "overlapping",
#'   lambda = lam.seq,
#'   group = pars.overlapping$groups,
#'   group_variable = pars.overlapping$groups_var,
#'   penalty_weights = pars.overlapping$group_weights,
#'   nfolds = 5,
#'   tol = 1e-4,
#'   maxit = 1e3,
#'   verbose = FALSE
#' )
#' 
#' plot(cv.overlapping)
#' plot(cv.overlapping, type = "solution-path")
#' 
#' @seealso \code{\link{sox}}, \code{\link{sox_cv}}.
#' @method plot sox_cv
#' @rdname plot.sox_cv
#' @export

plot.sox_cv <- function(x, type = "cv-curve", ...) {
  cvup <- x$cvup
  cvlo <- x$cvlo
  cvm <- x$cvm
  
  lam.seq <- x$lambdas
  lambda.min <- x$lambda.min
  lambda.1se <- x$lambda.1se
  
  if (type == "cv-curve") {
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
    axis(side = 3, at = log(lam.seq), labels = x$nzero, tick = FALSE)
  } else if (type == "solution-path") {
    estimates <- x$sox.fit$estimates
    
    matplot(x = log(lam.seq),
            y = t(estimates),
            xlab = expression(paste("Log(", lambda, ")")),
            ylab = "",
            type = "l")
    abline(h = 0, lty = 2)
    abline(v = log(lambda.min), lty = 3)
    abline(v = log(lambda.1se), lty = 3)
  }
}
