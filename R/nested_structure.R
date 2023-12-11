#' Automatically generate objects used to describe the structure of the nested group lasso penalty.
#'
#' @description
#' Automatically generate objects used to describe the structure of the nested group lasso penalty. The output is then used by \code{\link{sox}()} and \code{\link{sox_cv}()}.
#' 
#' @param group_list A list containing the indices of the group members.
#' 
#' @examples 
#' # p = 9 Variables:
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
#' # G = 12 Nested groups (misspecified, for the demonstration of the software only.)
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
#' str(pars.nested)
#'                 
#' @return A list of objects describing the group structure.
#'   \item{groups}{Required by \code{\link{sox}()} and \code{\link{sox_cv}()} to describe the relationship between the \eqn{G} \code{overlapping} groups. A \eqn{G * G} integer matrix whose \eqn{(i,j)} entry is \code{1} if and only if \eqn{i\neq j} and \eqn{g_i} is a child group (subset) of \eqn{g_j}, and is \code{0} otherwise.}
#'   \item{own_variables}{Required by \code{\link{sox}()} and \code{\link{sox_cv}()} to describe the relationship between the \eqn{G} \code{overlapping} groups and the \eqn{p} variables. The entries are the smallest variable indices in the groups (to achieve this, \code{group} is sorted. For any two groups \eqn{i} and \eqn{j}, if \eqn{i} is the parent group of \eqn{j}, then \eqn{i} is before \eqn{j} and vice versa, otherwise, the one with the smallest variable index is before the other.}
#'   \item{N_own_variables}{Required by \code{\link{sox}()} and \code{\link{sox_cv}()} to describe the relationship between the \eqn{G} \code{overlapping} groups and the \eqn{p} variables. An integer vector of length \eqn{G} indicating the number of variables that are in each group but not in any of its child groups.}
#'   \item{group_weights}{Required by \code{\link{sox}()} and \code{\link{sox_cv}()} to specify the group-specific penalty weights. The weight is generated in a way such that, the penalty weights of all the groups that contain a given variable sum to 1 for all variables.}

nested_structure <- function(group_list) {
  n_nodes <- length(group_list)
  
  group_list <- lapply(group_list, sort)
  
  reorder1 <- order(sapply(group_list, FUN = length), decreasing = T)
  
  group_list <- group_list[reorder1]
  
  reorder2 <- order(sapply(group_list, FUN = utils::head, n = 1))
  
  group_list <- group_list[reorder2]
  
  # pre-allocate own_variable, N_own_variable, group weights, groups matrix
  n_own_var <- integer(n_nodes)
  eta_g <- numeric(n_nodes)
  g <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  
  # own_variables
  own_var <- sapply(group_list, FUN = utils::head, n = 1)
  
  # N_own_variables
  for (i in seq(n_nodes)) {
    nodes_set <- group_list[-(1:i)]
    node_tmp <- group_list[[i]]
    nov_tmp <- 0L
    for (j in seq(length(node_tmp))) {
      if (!node_tmp[j] %in% unlist(nodes_set)) {
        nov_tmp <- nov_tmp + 1L
      }
    }
    n_own_var[i] <- nov_tmp
  }
  
  # groups
  for (i in 2:n_nodes) {
    for (j in (i-1):1) {
      if ( setequal(intersect(group_list[[i]], group_list[[j]]),  group_list[[i]]) ) {
        g[i, j] <- 1L
        break
      }
    }
  }
  
  # group weights: eta_g
  
  w <- pass_down <- numeric(n_nodes)
  pass_down[1] <- 1
  w[1] <- .Machine$double.eps
  
  sub_tree_old <- group_list
  
  for (i in 2:n_nodes) {
    node <- group_list[[i]]
    
    ancester <- which(g[i,] == 1L)
    
    sub_tree <- group_list[sapply(group_list, FUN = function(x) {setequal(x, intersect(x, node))})]
    sub_tree_ancester <- sub_tree
    sub_tree_ancester[[length(sub_tree_ancester) + 1]] <- group_list[[ancester]]
    
    ancester_height <- max(table(unlist(sub_tree_ancester)))
    node_height <- max(table(unlist(sub_tree)))
    pass_down[i] <- (node_height - 1) / (ancester_height - 1) * pass_down[ancester]
    
    if (ancester == 1) {
      w[i] <- 1 - pass_down[i]
    } else {
      w[i] <- w[ancester] / (1 - pass_down[ancester]) * pass_down[ancester] * (1 - pass_down[i])
    }
  }
  return(list(groups = g,
              own_variables = own_var,
              N_own_variables = n_own_var,
              group_weights = w))
}
