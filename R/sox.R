#' (Time-dependent) Cox model with structured variable selection
#'
#' @description
#' Fit a (time-dependent) Cox model with overlapping (including nested) group lasso penalty. The regularization path is computed at a grid of values for the regularization parameter lambda.
#' 
#' @param x Predictor matrix with dimension \eqn{nm * p}, where \eqn{n} is the number of subjects, \eqn{m} is the maximum observation time, and \eqn{p} is the number of predictors. See Details.
#' @param ID The ID of each subjects, each subject has one ID (multiple rows in \code{x} can share one \code{ID}).
#' @param time Represents the start of each time interval.
#' @param time2 Represents the stop of each time interval.
#' @param event Indicator of event. \code{event = 1} when event occurs and \code{event = 0} otherwise.
#' @param penalty Character string, indicating whether "\code{overlapping}" or "\code{nested}" group lasso penalty is imposed.
#' @param lambda Sequence of regularization coefficients \eqn{\lambda}'s.
#' @param group A \eqn{G * G} integer matrix required to describe the structure of the \code{overlapping} and \code{nested} groups. We recommend that the users generate it automatically using \code{\link{overlap_structure}()} and \code{\link{nested_structure}()}. See Examples and Details.
#' @param group_variable A \eqn{p * G} integer matrix required to describe the structure of the \code{overlapping} groups. We recommend that the users generate it automatically using \code{\link{overlap_structure}()}. See Examples and Details.
#' @param own_variable A non-decreasing integer vector of length \eqn{G} required to describe the structure of the \code{nested} groups. We recommend that the users generate it automatically using \code{\link{nested_structure}()}. See Examples and Details.
#' @param no_own_variable An integer vector of length \eqn{G} required to describe the structure of the \code{nested} groups. We recommend that the users generate it automatically using \code{\link{nested_structure}()}. See Examples and Details
#' @param penalty_weights Optional, vector of length \eqn{G} specifying the group-specific penalty weights. We recommend that the users generate it automatically using \code{\link{overlap_structure}()} or \code{\link{nested_structure}()}. If not specified, \eqn{\mathbf{1}_G} is used.
#' @param par_init Optional, vector of initial values of the optimization algorithm. Default initial value is zero for all \eqn{p} variables.
#' @param stepsize_init Initial value of the stepsize of the optimization algorithm. Default is 1.0.
#' @param stepsize_shrink Factor in \eqn{(0,1)} by which the stepsize shrinks in the backtracking linesearch. Default is 0.8.
#' @param tol Convergence criterion. Algorithm stops when the \eqn{l_2} norm of the difference between two consecutive updates is smaller than \code{tol}.
#' @param maxit Maximum number of iterations allowed.
#' @param verbose Logical, whether progress is printed.
#' 
#' @examples 
#' x <- as.matrix(sim[, c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B")])
#' lam.seq <- exp(seq(log(1e0), log(1e-3), length.out = 20))
#' 
#' # Variables:
#' ## 1: A1
#' ## 2: A2
#' ## 3: C1
#' ## 4: C2
#' ## 5: B
#' ## 6: A1B
#' ## 7: A2B
#' ## 8: C1B
#' ## 9: C2B
#' 
#' # Overlapping groups:
#' ## g1: A1, A2, A1B, A2B
#' ## g2: B, A1B, A2B, C1B, C2B
#' ## g3: A1B, A2B
#' ## g4: C1, C2, C1B, C2B
#' ## g5: C1B, C2B
#' 
#' overlapping.groups <- list(c(1, 2, 6, 7),
#'                            c(5, 6, 7, 8, 9),
#'                            c(6, 7),
#'                            c(3, 4, 8, 9),
#'                            c(8, 9))
#'                            
#' pars.overlapping <- overlap_structure(overlapping.groups)
#' 
#' fit.overlapping <- sox(
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
#'   tol = 1e-4,
#'   maxit = 1e3,
#'   verbose = FALSE
#' )
#' 
#' str(fit.overlapping)
#' 
#' # Nested groups (misspecified, for the demonstration of the software only.)
#' ## g1: A1, A2, C1, C2, B, A1B, A2B, C1B, C2B
#' ## g2: A1B, A2B, A1B, A2B
#' ## g3: C1, C2, C1B, C2B
#' ## g4: 1
#' ## g5: 2
#' ## ...
#' ## G12: 9
#' 
#' nested.groups <- list(1:9,
#'                       c(1, 2, 6, 7),
#'                       c(3, 4, 8, 9),
#'                       1, 2, 3, 4, 5, 6, 7, 8, 9)
#' 
#' pars.nested <- nested_structure(nested.groups)
#' 
#' fit.nested <- sox(
#'   x = x,
#'   ID = sim$Id,
#'   time = sim$Start,
#'   time2 = sim$Stop,
#'   event = sim$Event,
#'   penalty = "nested",
#'   lambda = lam.seq,
#'   group = pars.nested$groups,
#'   own_variable = pars.nested$own_variables,
#'   no_own_variable = pars.nested$N_own_variables,
#'   penalty_weights = pars.nested$group_weights,
#'   tol = 1e-4,
#'   maxit = 1e3,
#'   verbose = FALSE
#' )
#' 
#' str(fit.nested)
#' 
#' @details
#' The predictor matrix should be of dimension \eqn{nm * p}. Each row records the values of covariates for one subject at one time, for example, the values at the day from \code{time} (Start) to \code{time2} (Stop). An example dataset \code{\link{sim}} is provided. The dataset has the format produced by the \code{R} package \pkg{PermAlgo}. 
#' The specification of the arguments \code{group}, \code{group_variable}, \code{own_variable} and \code{no_own_variable} for the grouping structure can be found in \url{https://thoth.inrialpes.fr/people/mairal/spams/doc-R/html/doc_spams006.html#sec26} and \url{https://thoth.inrialpes.fr/people/mairal/spams/doc-R/html/doc_spams006.html#sec27}.
#'
#' In the Examples below, \eqn{p=9,G=5}, the group structure is: \deqn{g_1 = \{A_{1}, A_{2}, A_{1}B, A_{2}B\},} \deqn{g_2  = \{B, A_{1}B, A_{2}B, C_{1}B, C_{2}B\},} \deqn{g_3  = \{A_{1}B, A_{2}B\},} \deqn{g_4  = \{C_1, C_2, C_{1}B, C_{2}B\},} \deqn{g_5  = \{C_{1}B, C_{2}B\}.}
#' 
#' where \eqn{g_3} is a subset of \eqn{g_1} and \eqn{g_2}, and \eqn{g_5} is a subset of \eqn{g_2} and \eqn{g_4}.
#' @return A list with the following three elements.
#'   \item{lambdas}{The user-specified regularization coefficients \code{lambda} sorted in decreasing order.}
#'   \item{estimates}{A matrix, with each column corresponding to the coefficient estimates at each \eqn{\lambda} in \code{lambdas}.}
#'   \item{iterations}{A vector of number of iterations it takes to converge at each \eqn{\lambda} in \code{lambdas}.}
sox <- function(
    x,
    ID,
    time,
    time2,
    event,
    penalty,
    lambda,
    group,
    group_variable,
    own_variable,
    no_own_variable,
    penalty_weights,
    par_init,
    stepsize_init = 1,
    stepsize_shrink = 0.8,
    tol = 1e-5,
    maxit = 1000L,
    verbose = FALSE
) {
  
  p <- ncol(x)
  if (missing(par_init)) {
    par_init <- rep(0, p)
  }
  
  G <- nrow(group)
  if (missing(penalty_weights)) {
    penalty_weights <- rep(1, G)
  }
  
  n <- length(unique(ID))
  
  lambdas <- sort(lambda, decreasing = TRUE)
  
  if (missing(own_variable)) {
    if (penalty == "overlapping") {
      own_variable <- integer()
    } else {
      stop("'own_variable' must be provided for 'nested' group lasso.")
    }
  }
  if (missing(no_own_variable)) {
    if (penalty == "overlapping") {
      no_own_variable <- integer()
    } else {
      stop("'no_own_variable' must be provided for 'nested' group lasso.")
    }
  }
  
  if (missing(group_variable)) {
    if (penalty == "nested") {
      group_variable <- as.matrix(integer())
    } else {
      stop("'group_variable' must be provided for 'overlapping' group lasso.")
    }
  }
  
  if (penalty == "overlapping") {
    regul <- "graph"
  } else if (penalty == "nested") {
    regul <- "tree-l2"
    own_variable <- own_variable - 1L
  }
  
  fit <- sox_cpp(x = x,
                 start = time,
                 stop = time2,
                 event = event,
                 n_unique = n,
                 regul = regul,
                 lam = lambdas,
                 grp = group,
                 grpV = group_variable,
                 own_var = own_variable,
                 N_own_var = no_own_variable,
                 etaG = penalty_weights,
                 init = par_init,
                 l_ld = l_ld,
                 init_stepsize = stepsize_init,
                 ls_shrink = stepsize_shrink,
                 partol = tol,
                 maxit = maxit,
                 verbose = verbose)
  
  estimates <- fit$Estimates
  colnames(estimates) <- lambdas
  iterations <- fit$Iterations
  
  results <- list(lambdas = lambdas,
                  estimates = estimates,
                  iterations = iterations)
  
  if ( sum(iterations == maxit) > 0 ) {
    warning("Maximum iterations reached at some lambdas. Check ` $iterations`.")
  }
  
  class(results) <- "sox"
  return(results)
}