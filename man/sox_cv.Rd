% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sox_cv.R
\name{sox_cv}
\alias{sox_cv}
\title{cross-validation for \code{sox}}
\usage{
sox_cv(
  x,
  ID,
  time,
  time2,
  event,
  lambda,
  group,
  group_variable,
  penalty_weights,
  par_init,
  nfolds = 10,
  stepsize_init = 1,
  stepsize_shrink = 0.8,
  tol = 1e-05,
  maxit = 1000L,
  verbose = FALSE
)
}
\arguments{
\item{x}{Predictor matrix with dimension \eqn{nm * p}, where \eqn{n} is the number of subjects, \eqn{m} is the maximum observation time, and \eqn{p} is the number of predictors. See Details.}

\item{ID}{The ID of each subjects, each subject has one ID (many rows in \code{x} share one \code{ID}).}

\item{time}{Represents the start of each time interval.}

\item{time2}{Represents the stop of each time interval.}

\item{event}{Indicator of event. \code{event = 1} when event occurs and \code{event = 0} otherwise.}

\item{lambda}{Sequence of regularization coefficients \eqn{\lambda}'s.}

\item{group}{\eqn{G * G} matrix describing the relationship between the groups of variables, where \eqn{G} represents the number of groups. Denote the \eqn{i}-th group of variables by \eqn{g_i}. The \eqn{(i,j)} entry is \code{1} if and only if \eqn{i\neq j} and \eqn{g_i} is a child group (subset) of \eqn{g_j}, and is \code{0} otherwise. See Examples and Details.}

\item{group_variable}{\eqn{p * G} matrix describing the relationship between the groups and the variables. The \eqn{(i,j)} entry is \code{1} if and only if variable \eqn{i} is in group \eqn{g_j}, but not in any child group of \eqn{g_j}, and is \code{0} otherwise. See Examples and Details.}

\item{penalty_weights}{Optional, vector of length \eqn{G} specifying the group-specific penalty weights. If not specified, the default value is \eqn{\mathbf{1}_G}. Modify with caution.}

\item{par_init}{Optional, vector of initial values of the optimization algorithm. Default initial value is zero for all \eqn{p} variables.}

\item{nfolds}{Optional, the folds of cross-validation. Default is 10.}

\item{stepsize_init}{Initial value of the stepsize of the optimization algorithm. Default is 1.}

\item{stepsize_shrink}{Factor in \eqn{(0,1)} by which the stepsize shrinks in the backtracking linesearch. Default is 0.8.}

\item{tol}{Convergence criterion. Algorithm stops when the \eqn{l_2} norm of the difference between two consecutive updates is smaller than \code{tol}.}

\item{maxit}{Maximum number of iterations allowed.}

\item{verbose}{Logical, whether progress is printed.}
}
\value{
A list.
  \item{lambdas}{A vector of lambda used for each cross-validation.}
  \item{cvm}{The cv error averaged across all folds for each lambda.}
  \item{cvsd}{The standard error of the cv error for each lambda.}
  \item{cvup}{The cv error plus its standard error for each lambda.}
  \item{cvlo}{The cv error minus its standard error for each lambda.}
  \item{nzero}{The number of non-zero coefficients at each lambda.}
  \item{sox.fit}{A sox fit for the full data at all lambdas.}
  \item{lambda.min}{The lambda such that the \code{cvm} reach its minimum.}
  \item{lambda.1se}{The maximum of lambda such that the \code{cvm} is less than the minimum the \code{cvup} (the minmum of \code{cvm} plus its standard error).}
}
\description{
Conduct cross-validation (cv) for \code{sox}.
}
\details{
For each lambda, 10 folds cross-validation (by default) is performed. The cv error is defined as follows. Suppose we perform \eqn{K}-fold cross-validation, denote \eqn{\hat{\beta}^{-k}} by the estimate obtained from the rest of \eqn{K-1} folds (training set). The error of the \eqn{k}-th fold (test set) is defined as \eqn{2(P-Q)} divided by \eqn{R}, where \eqn{P} is the log partial likelihood evaluated at  \eqn{\hat{\beta}^{-k}} using the entire dataset, Q is the log partial likelihood evaluated at \eqn{\hat{\beta}^{-k}} using the training set, and R is the number of events in the test set. We do not use the negative log partial likelihood evaluated at \eqn{\hat{\beta}^{-k}} using the test set because the former definition can efficiently use the risk set, and thus it is more stable when the number of events in each test set is small (think of leave-one-out). The cv error is used in parameter tuning. To account for balance in outcomes among the randomly formed test set, we divide the deviance \eqn{2(P-Q)} by R. 
To get the estimated coefficients that has the minimum cv error, use \code{sox()$Estimates[sox()$Lambdas==sox_cv()$lambda.min]}. To apply the 1-se rule, use \code{sox()$Estimates[sox()$Lambdas==sox_cv()$lambda.1se]}.
}
\examples{
grp <- matrix(c(0, 0, 0, 0, 0,
                0, 0, 0, 0, 0,
                1, 1, 0, 0, 0,
                0, 0, 0, 0, 0,
                0, 1, 0, 1, 0),
              ncol = 5, byrow = TRUE)
grp.var <- matrix(c(1, 0, 0, 0, 0, #A1
                    1, 0, 0, 0, 0, #A2
                    0, 0, 0, 1, 0, #C1
                    0, 0, 0, 1, 0, #C2
                    0, 1, 0, 0, 0, #B
                    0, 0, 1, 0, 0, #A1B
                    0, 0, 1, 0, 0, #A2B
                    0, 0, 0, 0, 1, #C1B
                    0, 0, 0, 0, 1  #C2B
                   ), ncol = 5, byrow = TRUE)
eta_g <- rep(1, 5)
x <- as.matrix(sim[, c("A1","A2","C1","C2","B",
                       "A1B","A2B","C1B","C2B")])
lam.seq <- 10^seq(0, -2, by = -0.2)

cv <- sox_cv(x = x,
                ID = sim$Id,
                time = sim$Start,
                time2 = sim$Stop,
                event = sim$Event,
                lambda = lam.seq,
                group = grp,
                group_variable = grp.var,
                penalty_weights = eta_g,
                nfolds = 5,
                tol = 1e-4,
                maxit = 1e3,
                verbose = FALSE)
                
}
\seealso{
\code{\link{sox}}, \code{\link{plot_sox_cv}}.
}
