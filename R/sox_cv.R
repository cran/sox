#' cross-validation for \code{sox}
#'
#' @description
#' Conduct cross-validation (cv) for \code{sox}.
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
#' @param nfolds Optional, the folds of cross-validation. Default is 10.
#' @param foldid Optional, user-specified vector indicating the cross-validation fold in which each observation should be included. Values in this vector should range from 1 to \code{nfolds}. If left unspecified, \code{sox} will randomly assign observations to folds
#' @param stepsize_init Initial value of the stepsize of the optimization algorithm. Default is 1.
#' @param stepsize_shrink Factor in \eqn{(0,1)} by which the stepsize shrinks in the backtracking linesearch. Default is 0.8.
#' @param tol Convergence criterion. Algorithm stops when the \eqn{l_2} norm of the difference between two consecutive updates is smaller than \code{tol}.
#' @param maxit Maximum number of iterations allowed.
#' @param verbose Logical, whether progress is printed.
#' 
#' @details
#' For each lambda, 10 folds cross-validation (by default) is performed. The cv error is defined as follows. Suppose we perform \eqn{K}-fold cross-validation, denote \eqn{\hat{\beta}^{-k}} by the estimate obtained from the rest of \eqn{K-1} folds (training set). The error of the \eqn{k}-th fold (test set) is defined as \eqn{2(P-Q)} divided by \eqn{R}, where \eqn{P} is the log partial likelihood evaluated at  \eqn{\hat{\beta}^{-k}} using the entire dataset, Q is the log partial likelihood evaluated at \eqn{\hat{\beta}^{-k}} using the training set, and R is the number of events in the test set. We do not use the negative log partial likelihood evaluated at \eqn{\hat{\beta}^{-k}} using the test set because the former definition can efficiently use the risk set, and thus it is more stable when the number of events in each test set is small (think of leave-one-out). The cv error is used in parameter tuning. To account for balance in outcomes among the randomly formed test set, we divide the deviance \eqn{2(P-Q)} by R. 
#' To get the estimated coefficients that has the minimum cv error, use \code{sox_cv()$Estimates[, sox_cv$index["min",]]}. To apply the 1-se rule, use \code{sox_cv()$Estimates[, sox_cv$index["1se",]]}.
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
#' str(cv.overlapping)
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
#' cv.nested <- sox_cv(
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
#'   nfolds = 5,
#'   tol = 1e-4,
#'   maxit = 1e3,
#'   verbose = FALSE
#' )
#' 
#' str(cv.nested)
#'                 
#' @return A list.
#'   \item{lambdas}{A vector of lambda used for each cross-validation.}
#'   \item{cvm}{The cv error averaged across all folds for each lambda.}
#'   \item{cvsd}{The standard error of the cv error for each lambda.}
#'   \item{cvup}{The cv error plus its standard error for each lambda.}
#'   \item{cvlo}{The cv error minus its standard error for each lambda.}
#'   \item{nzero}{The number of non-zero coefficients at each lambda.}
#'   \item{sox.fit}{A fitted model for the full data at all lambdas of class "\code{sox}".}
#'   \item{lambda.min}{The lambda such that the \code{cvm} reach its minimum.}
#'   \item{lambda.1se}{The maximum of lambda such that the \code{cvm} is less than the minimum the \code{cvup} (the minmum of \code{cvm} plus its standard error).}
#'   \item{foldid}{The fold assignments used.}
#'   \item{index}{A one column matrix with the indices of \code{lambda.min} and \code{lambda.1se}.}
#'   \item{iterations}{A vector of number of iterations it takes to converge at each \eqn{\lambda} in \code{lambdas}.}
#' @seealso \code{\link{sox}}, \code{\link{plot.sox_cv}}.

sox_cv <- function(
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
    nfolds = 10,
    foldid = NULL,
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
  
  unique_ID <- unique(ID)
  n <- length(unique_ID)
  
  nlam <- length(lambda)
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
  
  ### sox fit
  
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
  
  sox.fit <- list(lambdas = lambdas,
                  estimates = estimates,
                  iterations = iterations)
  
  nzero <- colSums(apply(estimates, MARGIN = 2, "!=", 0))
  
  ### train-test split and pre-allocate result storage
  if (is.null(foldid)) {
    fold_id <- sample(rep(seq(nfolds), length = n))
  } else {
    fold_id <- foldid
  }
  cv_error <- matrix(NA, nrow = nfolds, ncol = nlam)
  iterations <- matrix(NA, nrow = nfolds, ncol = nlam)
  
  for ( k in seq(nfolds) ) {
    message("fold ", k)
    
    which_train <- unique_ID[fold_id != k]
    rows_train <- ID %in% which_train
    
    n_tmp <- length(which_train)
    
    train_x <- x[rows_train, ]
    train_time <- time[rows_train]
    train_time2 <- time2[rows_train]
    train_event <- event[rows_train]
    test_x <- x[!rows_train, ]
    test_time <- time[!rows_train]
    test_time2 <- time2[!rows_train]
    test_event <- event[!rows_train]
    
    fit <- sox_cpp(x = train_x,
                   start = train_time,
                   stop = train_time2,
                   event = train_event,
                   n_unique = n_tmp,
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
    
    deviance_full <- glmnet::coxnet.deviance(y = Surv(time,
                                                      time2,
                                                      event),
                                             x = x,
                                             beta = fit$Estimates)
    deviance_train <- glmnet::coxnet.deviance(y = Surv(train_time,
                                                       train_time2,
                                                       train_event),
                                              x = train_x,
                                              beta = fit$Estimates)
    
    cv_error[k, ] <- ( deviance_full - deviance_train )/sum(test_event)
    iterations[k, ] <- fit$Iterations
  }
  
  ### cv evaluation
  cvm <- apply(cv_error, MARGIN = 2, FUN = mean)
  cvsd <- apply(cv_error, MARGIN = 2, FUN = sd)/sqrt(nfolds)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  lam.min <- which.min(cvm)
  lam.1se <- sum(cvm[seq(lam.min)] > cvup[lam.min]) + 1
  lam.min.1se.id <- matrix(nrow = 2, ncol = 1)
  lam.min.1se.id[1, ] <- lam.min
  lam.min.1se.id[2, ] <- lam.1se
  rownames(lam.min.1se.id) <- c("min", "1se")
  colnames(lam.min.1se.id) <- "Lambda"
  
  results <- list(lambdas = lambdas,
                  cvm = cvm,
                  cvsd = cvsd,
                  cvup = cvup,
                  cvlo = cvlo,
                  nzero = nzero,
                  sox.fit = sox.fit,
                  lambda.min = lambdas[lam.min],
                  lambda.1se = lambdas[lam.1se],
                  foldid = fold_id,
                  index = lam.min.1se.id,
                  iterations = iterations)
  
  if ( sum(iterations == maxit) > 0 ) {
    warning("Maximum iterations reached at some lambdas. Check ` $iterations`.")
  }
  
  class(results$sox.fit) <- "sox"
  class(results) <- "sox_cv"
  return(results)
}

