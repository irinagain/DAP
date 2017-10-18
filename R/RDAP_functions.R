# Some modifications and wrappers for RDA based on sparsediscrim package from Ramey
# Purpose - take into account variables as well as lambda/gamma values, and have consistent folds between each

# This is a copy of update_hdrda function which for some reason is hidden
update_hdrda <- function (obj, lambda = 1, gamma = 0) 
{
  if (lambda == 0 && gamma == 0) {
    Gamma <- matrix(0, nrow = obj$q, ncol = obj$q)
  }
  else {
    Gamma <- obj$est[[1]]$alpha * lambda * obj$D_q + gamma
    Gamma_inv <- diag(Gamma^(-1))
  }
  for (k in seq_len(obj$num_groups)) {
    n_k <- obj$est[[k]]$n_k
    if (lambda == 0 && gamma == 0) {
      Q <- diag(n_k)
      W_k <- cov_mle(obj$est[[k]]$XU)
      W_inv <- try(solve_chol(W_k), silent = TRUE)
      if (inherits(W_inv, "try-error")) {
        W_k <- W_k + diag(0.001, nrow = nrow(W_k), ncol = ncol(W_k))
        W_inv <- solve_chol(W_k)
      }
    }
    else {
      XU_Gamma_inv <- obj$est[[k]]$XU %*% Gamma_inv
      Q <- diag(n_k) + obj$est[[k]]$alpha * (1 - lambda)/n_k * 
        tcrossprod(XU_Gamma_inv, obj$est[[k]]$XU)
      W_inv <- obj$est[[k]]$alpha * (1 - lambda)/n_k * 
        crossprod(XU_Gamma_inv, solve(Q, XU_Gamma_inv))
      W_inv <- Gamma_inv - W_inv
    }
    obj$est[[k]]$Gamma <- Gamma
    obj$est[[k]]$Q <- Q
    obj$est[[k]]$W_inv <- W_inv
  }
  obj$lambda <- lambda
  obj$gamma <- gamma
  obj
}

# This is a copy of predict.hdrda function which for some reason is hidden
predict.hdrda <- function (object, newdata, projected = FALSE, ...) 
{
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }
  newdata <- as.matrix(newdata)
  scores <- sapply(object$est, function(class_est) {
    if (object$lambda == 0 && object$gamma == 0) {
      log_det <- -log_determinant(class_est$W_inv)
    }
    else {
      log_det <- log_determinant(class_est$Q)
    }
    if (projected) {
      U1_x <- scale(newdata, center = class_est$xbar_U1, 
                    scale = FALSE)
      quad_forms <- diag(drop(tcrossprod(U1_x %*% class_est$W_inv, 
                                         U1_x)))
    }
    else {
      x_centered <- scale(newdata, center = class_est$xbar, 
                          scale = FALSE)
      U1_x <- crossprod(object$U1, t(x_centered))
      quad_forms <- apply(U1_x, 2, function(z) {
        quadform(class_est$W_inv, z)
      })
    }
    quad_forms + log_det - 2 * log(class_est$prior)
  })
  if (is.vector(scores)) {
    names(scores) <- object$groups
    min_scores <- which.min(scores)
    posterior <- exp(-(scores - min(scores)))
    posterior <- posterior/sum(posterior)
  }
  else {
    min_scores <- apply(scores, 1, which.min)
    posterior <- exp(-(scores - apply(scores, 1, min)))
    posterior <- posterior/rowSums(posterior)
  }
  class <- with(object, factor(groups[min_scores], levels = groups))
  list(class = class, scores = scores, posterior = posterior)
}

hdrda_cv_variables <- function (x, y, variables_list, num_folds = 5, num_lambda = 21, num_gamma = 8, 
          shrinkage_type = c("ridge", "convex"), verbose = FALSE, ...) 
{
  x <- as.matrix(x)
  y <- as.factor(y)
  n_var <- length(variables_list) #added
  shrinkage_type <- match.arg(shrinkage_type)
  cv_folds <- cv_partition(y = y, num_folds = num_folds)
  seq_lambda <- seq(0, 1, length = num_lambda)
  if (shrinkage_type == "ridge") {
    seq_gamma <- c(0, 10^seq.int(-2, num_gamma - 4))
  } else {
    seq_gamma <- seq(0, 1, length = num_gamma)
  }
  tuning_grid <- expand.grid(lambda = seq_lambda, gamma = seq_gamma) 
  tuning_grid <- dplyr::arrange(tuning_grid, lambda, gamma)
  errors_var <- matrix(0, nrow(tuning_grid), n_var) #added to store errors over variables
  for (i in 1:n_var){ #added loop over variables
    cv_errors <- sapply(seq_along(cv_folds), function(fold_i) {
      if (verbose) {
        cat("CV Fold: ", fold_i, " of ", num_folds, "\n")
      }
      fold <- cv_folds[[fold_i]]
      train_x <- x[fold$training, variables_list[[i]], drop = F]
      train_y <- y[fold$training]
      test_x <- x[fold$test, variables_list[[i]], drop = F]
      test_y <- y[fold$test]
      hdrda_out <- hdrda(x = train_x, y = train_y, lambda = 1, 
                       gamma = 0, shrinkage_type = shrinkage_type)
      test_x <- test_x %*% hdrda_out$U1
      fold_errors <- mapply(function(lambda, gamma) {
      errors <- try({
        hdrda_updated <- update_hdrda(hdrda_out, lambda, 
                                      gamma)
        sum(predict.hdrda(hdrda_updated, test_x, projected = TRUE)$class != 
              test_y)
      }, silent = TRUE)
      errors
    }, tuning_grid$lambda, tuning_grid$gamma)
    fold_errors
    })
    errors_var[,i] <- rowSums(cv_errors)/nrow(x) # added to store errors over variables
  }
  cv_summary <- cbind(tuning_grid, errors_var) #modified
  optimal <- which(errors_var==min(errors_var), arr.ind = T) #modified
  optimal <- optimal[1,] # always choose the first
  
  lambda <- tuning_grid$lambda[optimal[1]] #modified
  gamma <- tuning_grid$gamma[optimal[1]] #modified
  variables <- variables_list[[optimal[2]]] #added
  hdrda_out <- hdrda(x = x[,variables, drop = F], y = y, lambda = lambda, gamma = gamma, 
                     shrinkage_type = shrinkage_type)
  hdrda_out$lambda <- lambda
  hdrda_out$gamma <- gamma
  hdrda_out$variables <- variables #added
  hdrda_out$cv_summary <- cv_summary
  class(hdrda_out) <- c("hdrda_cv", "hdrda")
  hdrda_out
}

quadform <- function (A, x){
  drop(crossprod(x, A %*% x))
}

# HDRDA as in Ramey which is actually similar to Friedman's RDA
apply_hdrda_full <- function(xtrain, ytrain, xtest, ytest){
  outcv <- hdrda_cv(x = xtrain, y = ytrain, shrinkage_type="convex", num_folds = 5)
  return(list(error = mean(predict(outcv,xtest)$class != ytest), features = ncol(xtrain), features_id = c(1:ncol(xtrain))))
}

apply_RDAP <- function(xtrain, ytrain, xtest, ytest, lambda_seq = NULL, n_lambda = 50,  maxmin_ratio = 0.1, eps = 1e-4, m_max = 10000){
  Xmean <- colMeans(xtrain)
  xtrain <- xtrain - matrix(Xmean, nrow(xtrain), ncol(xtrain), byrow = T)
  xtest <- xtest - matrix(Xmean, nrow(xtest), ncol(xtest), byrow = T)
  out_s <- standardizeData(xtrain, ytrain, center = F)
  
  ## generate lambda sequence
  l_max <- max(sqrt(colMeans(out_s$X1)^2 + colMeans(out_s$X2)^2))
  if (!is.null(lambda_seq)) {
    nl= length(lambda_seq)
    if (nl < 1){
      warning(paste("There is no qualified lambda value. New values will be generated automatically. n_lambda will be set as.", n_lambda,sep = " "))
      lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
    }else{
      n_lambda = nl
    }  
  }else {
    lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
  }
  lambda_seq = sort(lambda_seq, decreasing = TRUE)
  
  # Solve the whole problem on the sequence of lamdba values
  fit = solve_DAP_seq(X1 = out_s$X1, X2 = out_s$X2, lambda_seq = lambda_seq, eps = eps, m_max = m_max)
  
  # Identify the variables list
  # Index of lambda values corresponding to non-zero solution
  nonzero = which(fit$nfeature_vec > 1)
  variables_list = list()
  iter = 1
  for (i in nonzero){
    if (iter == 1){
      variables_list[[iter]] = which((abs(fit$V1_mat[,i]) + abs(fit$V2_mat[,i]))>0)
      var_old = variables_list[[iter]]
      iter = iter + 1
    }else{
      var_new = which((abs(fit$V1_mat[,i]) + abs(fit$V2_mat[,i]))>0)
      # Remove repetitions
      if (length(var_new) != length(var_old)){
        variables_list[[iter]] = var_new
        var_old = var_new
        iter = iter + 1
      }
    }
  }
  
  # Apply cv
  outcv <- hdrda_cv_variables(x = xtrain, y = ytrain, variables_list = variables_list, shrinkage_type="convex")
  return(list(error = mean(predict(outcv,xtest[,outcv$variables])$class != ytest), features = length(outcv$variables), features_id = outcv$variables))
}