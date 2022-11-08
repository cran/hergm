
## This file contains helper functions

norm2 <- function(vec) { 
  val <- norm(as.matrix(vec), type = "2")
  return(val)
}

norm1 <- function(vec) { 
  val <- norm(as.matrix(vec), type = "1")
  return(val)
}


rep_row <- function(x, n) {
       matrix(rep(x, each = n), nrow = n)
}


exp_fun <- function(x) {
  val <- exp(x)
  if (length(val) > 1) {
    for (i in 1:length(val)) {
      if (val[i] > .Machine$double.xmax) {
        val[i] <- .Machine$double.xmax
      }
    }
  } else {
   if (val > .Machine$double.xmax) {
     val <- .Machine$double.xmax
   }
 }
 return(val)
}


log_fun <- function(x) {
  val <- log(x)
  if (val > .Machine$double.xmax) {
    val <- .Machine$double.xmax
  }
  if (val < .Machine$double.xmin) {
    val <- .Machine$double.xmin
  }
  return(val)
}


outer_fun <- function(v) {
  outer(v, v)
}


# Sort and order the neighborhoods by size
sort_neighborhoods <- function(net_list, mod_names) {

  # Get the sizes of each neighborhood
  n_ <- numeric(length(net_list))
  for (i in 1:length(net_list)) {
    n_[i] <- net_list[[i]]$gal$n
  }

  sort_order <- order(n_)
  net_list_order <- rep(list(NULL), length(net_list))
  mod_names_order <- rep(list(NULL), length(mod_names))
  for (i in 1:length(net_list)) {
    net_list_order[[i]] <- net_list[[sort_order[i]]]
    mod_names_order[[i]] <- mod_names[[sort_order[i]]]
  }
  return(list(net_list = net_list_order, mod_names = mod_names_order, sort_order = sort_order))
}

dim_fun <- function(n_obj, n_groups, eta_len) {
  sizes <- split_neighborhoods(n_obj, n_groups)
  dims <- matrix(0, nrow = n_obj, ncol = eta_len)
  for (i in 1:n_groups) {
    num_ <- sizes[i]
    dim_1 <- 1 + sum(sizes[seq_spec(i, adjust = -1)])
    dim_2 <- sum(sizes[1:i])
    dims[dim_1:dim_2, ] <-
          rep_row(seq(1 + (i - 1) * eta_len, eta_len * i), sizes[i])
  }
  return(list(dims = dims, sizes = sizes))
}

split_neighborhoods <- function(n_obj, n_groups) {
  base_ <- n_obj %/% n_groups
  left_ <- n_obj %% n_groups
  sizes <- rep(base_, n_groups) + c(rep(1, left_), rep(0, n_groups - left_))
  return(sizes)
}

seq_spec <- function(i, adjust = 0) {
  if (i == 1) {
    return(numeric(0))
  } else {
    return(1:(i + adjust))
  }
}

make_eta_fun <- function(num_group, parameterization) {
  if (parameterization == "multi_group") {
    eta_fun <- function(eta) {
      num_ <- 1
      len_ <- length(eta) / num_
      eta_base <- eta[1:len_]
      eta_val <- eta_base
      for (i in 2:num_) {
        dim_1 <- 1 + len_ * (i - 1)
        dim_2 <- len_ * i
        cur_val <- eta_base + eta[dim_1:dim_2]
        eta_val <- c(eta_val, cur_val)
      }
      return(eta_val)
    }
    body(eta_fun)[[2]] <- substitute(num_ <- num_group,
                                     list(num_group = num_group))

  } else if (parameterization == "size") {
    eta_fun <- function(eta) {
      return(eta)
    }
  }
  return(eta_fun)
}


make_eta_grad <- function(num_group, parameterization) {
  if (parameterization == "multi_group") {
    eta_grad <- function(eta) {
      num_ <- 1
      len_ <- length(eta) / num_
      eta_grad_val <- diag(len_)
      for (i in 2:num_) {
        eta_grad_val <- as.matrix(bdiag(eta_grad_val, diag(len_)))
      }
      eta_grad_val[ , 1:len_] <- rbind(t(matrix(rep(diag(len_), num_group),
                                         nrow = len_,
                                         ncol = num_group * len_)))
      return(eta_grad_val)
    }
    body(eta_grad)[[2]] <- substitute(num_ <- num_group,
                                      list(num_group = num_group))

  } else if (parameterization == "size") {
    eta_grad <- function(eta) {
      return(diag(length(eta)))
    }
  }
  return(eta_grad)
}


assign_labels <- function(K, sizes) {
  labels <- numeric(K)
  size_ <- c(0, sizes)
  for (i in 1:K) {
    labels[i] <- max(which(i > cumsum(size_)))
  }
  return(labels)
}


make_return_obj <- function(obj, labels, sort_order) {
  n_ <- length(unique(labels))
  return_list <- rep(list(NULL), n_)
  len_ <- length(obj$est$eta) / n_
  # names(return_list) <- Rprintf("group%i", 1:n_)
  grad <- obj$est$eta_grad(obj$est$eta)
  info_mat <- t(solve(grad)) %*% obj$est$info_mat %*% solve(grad)
  se_vec <- sqrt(diag(solve(info_mat)))
  for (i in 1:n_) {
    return_list[[i]] <- list(labels = NULL, estimates = NULL, se = NULL)
    return_list[[i]]$labels <- sort(sort_order[labels == i])

    dim_1 <- 1 + len_ * (i - 1)
    dim_2 <- len_ * i
    return_list[[i]]$estimates <- obj$est$eta_fun(obj$est$eta)[dim_1:dim_2]
    return_list[[i]]$se <- se_vec[dim_1:dim_2]
  }
  return(return_list)
}


check_extensions <- function(mod_names) { 
  L <- length(mod_names) 
  for (i in 1:L) { 
    mod_names[[i]] <- strsplit(as.character(mod_names[[i]]), "_ijk")
    mod_names[[i]] <- strsplit(as.character(mod_names[[i]]), "_ij") 
  }
  return(mod_names) 
}
