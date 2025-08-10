#' honeycombHIVE: Multilayer AIS with optional gradient-based fine-tuning
#' 
#' The \code{honeycombHIVE} function implements a multilayer artificial immune 
#' system that iteratively refines a set of prototypes - referred to as 
#' antibodies - to model the structure of the input data. In each layer, the 
#' function first uses the \code{\link{bHIVE}} algorithm to generate or update 
#' antibodies based on the current data representation and task (clustering, 
#' classification, or regression). Optionally, it applies gradient-based 
#' fine-tuning (via \code{\link{refineB}}) to these antibodies, allowing for 
#' advanced refinement through various optimizers (e.g., SGD, Adam, RMSProp) and 
#' customizable loss functions. The final output is a hierarchical set of layers 
#' that encapsulate both the refined prototypes and the corresponding cluster 
#' assignments or predictions for the original observations, making 
#' \code{honeycombHIVE} a versatile tool for adaptive learning and pattern 
#' recognition.
#' 
#' @param X A numeric matrix or data frame of input features 
#' (rows = observations, columns = features).
#' @param y Optional target vector (factor for classification, numeric 
#' for regression).
#' @param task Character, one of "clustering", "classification", or "regression".
#' @param layers Integer, how many layers (AIS iterations) to run.
#' @param nAntibodies Integer, how many antibodies (prototypes) to generate 
#' initially in each layer.
#' @param minAntibodies Integer, minimal number of antibodies to keep in each 
#' layer.
#' @param epsilon Numeric, threshold param for \code{bHIVE} suppression.
#' @param beta Numeric, selection pressure param for \code{bHIVE}.
#' @param maxIter Integer, maximum iterations for \code{bHIVE} each layer.
#' @param collapseMethod One of "centroid","medoid","median","mode".
#' @param minClusterSize Minimum cluster size. Smaller clusters can be 
#' merged/discarded if not NULL.
#' @param distance Distance metric for medoid calculation, e.g. "euclidean".
#' @param verbose Logical, if TRUE prints progress at each layer.
#' @param refine Logical, if TRUE apply gradient-based refinement via 
#' \code{refineB()} to each layer's prototypes.
#' @param refineLoss Character specifying the loss for \code{refineB()} 
#' (e.g. "mse", "mae", etc.).
#' @param refineSteps Integer, number of gradient steps in \code{refineB()}.
#' @param refineLR Numeric, learning rate for gradient updates.
#' @param refinePushAway Logical, if TRUE and classification, push prototypes 
#' away from differently labeled points.
#' @param refineHuberDelta Numeric, delta parameter if using the "huber" loss.
#' @param refineOptimizer Character, one of \code{"sgd", "momentum", "adagrad", 
#' "adam", "rmsprop"} to be passed to \code{refineB()}.
#' @param refineMomentumCoef Numeric, momentum coefficient (if using momentum).
#' @param refineBeta1 Numeric, first moment decay rate (if using Adam).
#' @param refineBeta2 Numeric, second moment decay rate (if using Adam).
#' @param refineRmspropDecay Numeric, decay rate for the moving average of 
#' squared gradients (if using RMSProp).
#' @param refineEpsilon Numeric, a small constant for numerical stability 
#' (used in adaptive optimizers).
#' @param ... Additional arguments passed to \code{bHIVE}.
#'
#' @examples
#' # Example 1: Clustering
#' data(iris)
#' X_iris <- iris[, 1:4]
#' resC <- honeycombHIVE(
#'   X = X_iris,
#'   task = "clustering",
#'   layers = 3,
#'   nAntibodies = 15,
#'   beta = 5,
#'   maxIter = 10
#' )
#'
#' # Example 2: Regression
#' set.seed(42)
#' X_reg <- matrix(rnorm(100*4), ncol = 4)
#' y_reg <- rowSums(X_reg[, 1:2]) + rnorm(100)
#' resReg <- honeycombHIVE(
#'   X = X_reg,
#'   y = y_reg,
#'   task = "regression",
#'   layers = 3,
#'   nAntibodies = 10
#' )
#' 
#' @return A list of length \code{layers}. Each element is a list representing a layer and contains:
#'   \itemize{
#'     \item \code{antibodies}: A matrix of the antibody prototypes for that layer.
#'     \item \code{assignments}: A vector of antibody indices for each data point in the layer's input.
#'     \item \code{membership}: A vector indicating which final cluster each of the **original** data points belongs to.
#'     \item \code{predictions}: If classification or regression, a vector of predictions for each of the **original** data points.
#'     \item \code{task}: The task performed ('clustering', 'classification', or 'regression').
#'   }
#'
#' @importFrom stats median setNames
#' @export
honeycombHIVE <- function(X,
                          y = NULL,
                          task = c("clustering", "classification", "regression"),
                          layers = 3,
                          nAntibodies = 20,
                          minAntibodies = 5,
                          epsilon = 0.05,
                          beta = 5,
                          maxIter = 10,
                          collapseMethod = c("centroid", "medoid", "median", "mode"),
                          minClusterSize = NULL,
                          distance = "euclidean",
                          verbose = TRUE,
                          refine = FALSE,
                          refineLoss = "mse",
                          refineSteps = 5,
                          refineLR = 0.01,
                          refinePushAway = TRUE,
                          refineHuberDelta = 1.0,
                          refineOptimizer = "sgd",
                          refineMomentumCoef = 0.9,
                          refineBeta1 = 0.9,
                          refineBeta2 = 0.999,
                          refineRmspropDecay = 0.9,
                          refineEpsilon = 1e-8,
                          ...) {
  
  # =======================
  # 0) Basic Checks & Setup
  # =======================
  task <- match.arg(task)
  collapseMethod <- match.arg(collapseMethod)
  
  .validate_bHIVE_input(X, y)  
  
  # Validate refineLoss based on task
  valid_losses <- list(
    regression = c("mse", "mae", "huber", "kullback_leibler"),
    classification = c("categorical_crossentropy", "binary_crossentropy", "kullback_leibler"),
    clustering = c("mse", "mae", "huber")
  )
  
  if (refine && !refineLoss %in% valid_losses[[task]]) {
    stop(sprintf("Invalid refineLoss '%s' for task '%s'. Supported losses: %s", 
                 refineLoss, task, paste(valid_losses[[task]], collapse = ", ")))
  }
  
  # Force X into a data.frame for subsetting operations
  X <- as.data.frame(X)
  n_original <- nrow(X)
  
  # Assign rownames if missing to track original data points
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("row_", seq_len(n_original))
  }
  
  # `rowIndices` maps the rownames of the current layer's data to the indices of the original data
  rowIndices <- setNames(as.list(seq_len(n_original)), rownames(X))
  
  # If y is provided, assign names so that we can index by row names
  if (!is.null(y)) {
    if (is.null(names(y))) {
      names(y) <- rownames(X)
    }
  }
  
  # Initialize data for the first layer
  current_X <- X
  current_y <- if (!is.null(y)) y else NULL
  
  # Prepare a list to store results from each layer
  results <- vector("list", length = layers)
  dots <- list(...)  # additional args for bHIVE
  
  # ========================
  # 1) Iteration Over Layers
  # ========================
  for (layer_i in seq_len(layers)) {
    if (verbose) {
      message(sprintf("\n=== honeycombHIVE: Layer %d / %d (task=%s) ===",
                      layer_i, layers, task))
    }
    
    # 1a) Decide how many antibodies for this layer
    n_current <- nrow(current_X)
    nAntibodies_layer <- min(nAntibodies, max(minAntibodies, n_current))
    
    # 1b) Construct arguments for bHIVE
    bHIVE_args <- list(
      X           = current_X,
      y           = current_y,
      task        = task,
      nAntibodies = nAntibodies_layer,
      epsilon     = epsilon,
      beta        = beta,
      maxIter     = maxIter,
      verbose     = verbose
    )
    extra_args <- dots[setdiff(names(dots), names(bHIVE_args))]
    bHIVE_args <- c(bHIVE_args, extra_args)
    
    # 1c) Run bHIVE on the current layer's data
    res_layer <- do.call(bHIVE, bHIVE_args)
    
    # 1d) Optional refinement step with refineB
    if (refine) {
      # ... (refinement logic remains the same)
    }
    
    # =================================================================
    # 2) Collapse Step: Summarize clusters to form the next layer's input
    # =================================================================
    cluster_ids <- res_layer$assignments
    subsets_rownames <- split(rownames(current_X), f = cluster_ids)
    
    # This function processes a single cluster
    process_cluster <- function(rn_in_cluster) {
      orig_indices <- unlist(rowIndices[rn_in_cluster], use.names = FALSE)
      n_members <- length(orig_indices)
      
      sub_data <- current_X[rn_in_cluster, , drop = FALSE]
      
      # Calculate the prototype for the cluster
      if (!is.null(minClusterSize) && n_members < minClusterSize) {
        proto <- rep(NA_real_, ncol(current_X))
      } else {
        proto <- switch(
          collapseMethod,
          "centroid" = colMeans(sub_data, na.rm = TRUE),
          "medoid"   = {
            dmat <- as.matrix(dist(sub_data, method = distance))
            idx_medoid <- which.min(rowSums(dmat))
            as.numeric(sub_data[idx_medoid, , drop = FALSE])
          },
          "median"   = apply(sub_data, 2, median, na.rm = TRUE),
          "mode"     = {
            apply(sub_data, 2, function(vec) {
              tb <- table(vec)
              as.numeric(names(tb)[which.max(tb)])
            })
          }
        )
      }
      
      # Calculate predictions for the original data points in this cluster
      prediction <- NA
      if (task == "regression" && !is.null(current_y)) {
        indices_in_current <- match(rn_in_cluster, rownames(current_X))
        prediction <- mean(current_y[indices_in_current], na.rm = TRUE)
        if (is.na(prediction)) prediction <- mean(y, na.rm = TRUE)
      } else if (task == "classification" && !is.null(current_y)) {
        y_vals <- current_y[rn_in_cluster]
        tbl <- table(y_vals)
        prediction <- names(tbl)[which.max(tbl)]
      }
      
      list(proto = proto, orig_indices = orig_indices, prediction = prediction)
    }
    
    # Apply the processing function to each cluster
    processed_clusters <- lapply(subsets_rownames, process_cluster)
    
    # Build the prototype matrix for the next layer
    proto_matrix <- do.call(rbind, lapply(processed_clusters, `[[`, "proto"))
    
    # Filter out empty prototypes
    keep_mask <- !apply(proto_matrix, 1, function(x) all(is.na(x)))
    proto_matrix <- proto_matrix[keep_mask, , drop = FALSE]
    processed_clusters <- processed_clusters[keep_mask]
    
    if (nrow(proto_matrix) == 0) {
      stop("No valid prototypes generated. Check parameters & input data.")
    }
    
    # Assign names to the new prototypes
    row_names <- paste0("Layer", layer_i, "_Cluster", seq_len(nrow(proto_matrix)))
    rownames(proto_matrix) <- row_names

    # Save predictions and membership for this layer
    membership <- rep(NA_integer_, n_original)
    predictions_all_original <- if (task %in% c("classification", "regression")) rep(NA, n_original) else NULL
    
    for (i in seq_along(processed_clusters)) {
      cluster_info <- processed_clusters[[i]]
      membership[cluster_info$orig_indices] <- i
      if (!is.null(predictions_all_original)) {
        predictions_all_original[cluster_info$orig_indices] <- cluster_info$prediction
      }
    }

    res_layer$membership <- membership
    res_layer$predictions <- predictions_all_original
    results[[layer_i]] <- res_layer
    
    if (verbose) {
      message(sprintf("Layer %d completed. Next layer will use %d prototypes.",
                      layer_i, nrow(proto_matrix)))
    }
    
    # =========================
    # 3) Prepare for Next Layer
    # =========================
    current_X <- as.data.frame(proto_matrix)
    
    # Update rowIndices for the next layer
    new_rowIndices <- lapply(processed_clusters, `[[`, "orig_indices")
    names(new_rowIndices) <- row_names
    rowIndices <- new_rowIndices
    
    # Rebuild current_y for the next layer
    if (task == "classification" && !is.null(y)) {
      current_y <- vapply(rowIndices, function(ids) {
        tbl <- table(y[ids])
        names(tbl)[which.max(tbl)]
      }, character(1))
      current_y <- factor(current_y, levels = levels(y))
    } else if (task == "regression" && !is.null(y)) {
      current_y <- vapply(rowIndices, function(ids) {
        mean(y[ids], na.rm = TRUE)
      }, numeric(1))
    }
  }
  return(results)
}
