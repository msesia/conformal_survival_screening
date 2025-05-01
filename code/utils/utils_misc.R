get_pretty_quantiles <- function(x, n = 10) {
  if (!is.numeric(x) || length(x) == 0) stop("x must be a non-empty numeric vector.")
  x <- x[!is.na(x) & is.finite(x) & x > 0]  # Filter to positive real values
  
  # Compute quantiles
  probs <- seq(0, 1, length.out = n)
  q_raw <- quantile(x, probs = probs, na.rm = TRUE)
  
  # Smart rounding for readability
  round_to_nice <- function(vals) {
    rng <- range(vals)
    spread <- diff(rng)
    
    if (spread > 100) {
      digits <- 0
    } else if (spread > 10) {
      digits <- 1
    } else if (spread > 1) {
      digits <- 2
    } else {
      digits <- 3
    }
    round(vals, digits)
  }
  
  q_rounded <- round_to_nice(q_raw)
  q_unique <- unique(q_rounded)
  
  # If we lost too many points due to rounding, fallback to less aggressive rounding
  if (length(q_unique) < n) {
    q_unique <- unique(round(q_raw, 4))
  }
  
  return(sort(q_unique))
}

evaluate_classification <- function(selections, oracle = NULL, n_test = 1000) {
  # If no oracle provided, use selections$oracle
  if (is.null(oracle)) {
    oracle <- selections$oracle
  }
  
  # Define the full set of indices (test set)
  all_indices <- 1:n_test
  
  # True positive (high-risk) and true negative sets
  true_positives <- oracle
  true_negatives <- setdiff(all_indices, oracle)
  
  # List all methods to evaluate
  methods <- names(selections)
  
  # Evaluate each method
  results <- lapply(methods, function(method) {
    selected <- selections[[method]]
    
    TP <- sum(selected %in% true_positives)
    FP <- sum(selected %in% true_negatives)
    FN <- sum(!(true_positives %in% selected))
    TN <- sum(!(true_negatives %in% selected))
    
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    f1_score  <- ifelse(!is.na(precision) && !is.na(recall) && (precision + recall) > 0,
                        2 * (precision * recall) / (precision + recall), NA)
    accuracy  <- (TP + TN) / (TP + FP + FN + TN)
    
    return(data.frame(
      Method = method,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      Precision = precision,
      Recall = recall,
      F1_score = f1_score,
      Accuracy = accuracy
    ))
  })
  
  # Bind all results into a nice data frame
  do.call(rbind, results)
}

