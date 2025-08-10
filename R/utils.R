utils::globalVariables(c("PC1", "PC2", "Group", "Feature", "Prototype", "Layer"))

.validate_bHIVE_input <- function(X, 
                                  y = NULL) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("Input X must be a matrix or data frame.")
  }
  if (!is.null(y)) {
    if (!is.factor(y) & !is.numeric(y)) {
      stop("y must be a factor (classification) or numeric vector (regression).")
    }
    if (nrow(X) != length(y)) {
      stop("X and y must have the same number of rows.")
    }
  }
  if (anyNA(X)) {
    stop("Input X contains missing values. Please handle missing values before running bHIVE.")
  }
  invisible(TRUE)
}

# Determines how to move a single prototype 'ab_vec' in feature space.
# This function calculates a heuristic gradient for a given loss function.
# NOTE: The gradient calculations, especially for classification, are simplified
# heuristics (push/pull) rather than exact mathematical gradients. This is due

# to the nature of the `refineB` function, which operates on assignments
# post-hoc without direct access to model probabilities.
.update_prototype <- function(ab_vec, 
                              x_i, 
                              same_label = NA, 
                              task = c("clustering","classification","regression"),
                              loss = c("categorical_crossentropy", "binary_crossentropy", 
                                       "kullback_leibler", "cosine",
                                       "mse","mae","poisson","huber"),
                              push_away = TRUE,
                              huber_delta = 1.0) {
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  # The gradient must be a numeric vector of the same length as the antibody
  d <- length(ab_vec)
  zero_vec <- numeric(d)          # A zero vector for quick returns
  pull_vec <- x_i - ab_vec        # Heuristic gradient to pull ab_vec toward x_i
  push_vec <- -pull_vec           # Heuristic gradient to push ab_vec away from x_i
  
  # =================
  # 1) CLUSTERING
  # =================
  if (task == "clustering") {
    # For clustering, the goal is to move the prototype towards the data point.
    return(pull_vec)
  }
  
  # =================
  # 2) CLASSIFICATION
  # =================
  if (task == "classification") {
    # If the antibody has no assigned label, do not move it.
    if (is.na(same_label)) {
      return(zero_vec)
    }
    is_same <- isTRUE(same_label)
    
    # The classification losses use a simplified push/pull heuristic.
    # (A) If the data point has the same label as the antibody, PULL the antibody towards the data point.
    # (B) If the labels are different and `push_away` is TRUE, PUSH the antibody away.
    # (C) Otherwise, do nothing.
    if (loss %in% c("categorical_crossentropy","binary_crossentropy","kullback_leibler","cosine","mse")) {
      # For these losses, the heuristic is a simple vector difference.
      if (is_same) {
        return(pull_vec)
      } else if (push_away) {
        return(push_vec)
      } else {
        return(zero_vec)
      }
    }
    
    else if (loss == "mae") {
      # For MAE, the heuristic is based on the sign of the vector difference.
      if (is_same) {
        return(sign(pull_vec))
      } else if (push_away) {
        return(sign(push_vec))
      } else {
        return(zero_vec)
      }
    }
    
    # For other losses not suitable for classification, do nothing.
    return(zero_vec)
  }
  
  # =================
  # 3) REGRESSION
  # =================
  if (task == "regression") {
    # The residual in the feature space.
    residual <- pull_vec  # i.e., x_i - ab_vec
    
    if (loss == "mse") {
      # The gradient of 0.5 * ||residual||^2 is the residual itself.
      return(residual)
    }
    else if (loss == "mae") {
      # The gradient of ||residual||_1 is the sign of the residual.
      return(sign(residual))
    }
    else if (loss == "poisson") {
      # This is a heuristic for Poisson loss. A true Poisson loss would
      # require a rate parameter, which is not available here. We fall back
      # to the MSE-like gradient.
      return(residual)
    }
    else if (loss == "huber") {
      # Huber loss combines L2 and L1 norms.
      dist_val <- sqrt(sum(residual^2))
      if (dist_val <= huber_delta) {
        # In the L2 region, the gradient is the residual.
        return(residual)
      } else {
        # In the L1 region, the gradient is scaled.
        return(huber_delta * sign(residual))
      }
    }
    # Fallback for other losses: use the MSE-like gradient.
    return(residual)
  }
  
  # Should not be reached.
  return(zero_vec)
}
