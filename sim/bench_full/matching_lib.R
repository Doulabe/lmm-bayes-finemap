# ==============================================================================
# matching_lib.R
# One-to-one maximum bipartite matching metrics for LD-aware scoring.
# A selected variant s can be matched to a causal variant c iff |r(s,c)| >= tau.
# TP = size of a maximum matching; FP = K_hat - TP; FN = K_star - TP.
# Kuhn augmenting-path algorithm (right side = causals, small).
# ==============================================================================

max_bipartite_matching <- function(A) {
  # A: logical matrix, rows = selected, cols = causal
  nL <- nrow(A); nR <- ncol(A)
  if (nL == 0L || nR == 0L) return(0L)
  matchL <- integer(nL)                       # 0 = unmatched, else causal index
  total <- 0L
  for (r in seq_len(nR)) {
    vis <- logical(nL)
    aug <- function(rr) {
      for (s in which(A[, rr])) {
        if (!vis[s]) {
          vis[s] <<- TRUE
          if (matchL[s] == 0L || aug(matchL[s])) {
            matchL[s] <<- rr
            return(TRUE)
          }
        }
      }
      FALSE
    }
    if (aug(r)) total <- total + 1L
  }
  total
}

matching_metrics <- function(sel, causal, R, tau = 0.5) {
  K_hat <- length(sel); K_star <- length(causal)
  if (K_hat == 0L)
    return(data.frame(K_hat = 0L, TP = 0L, recall = 0, precision = 0, f1 = 0))
  A <- abs(R[sel, causal, drop = FALSE]) >= tau
  TP <- max_bipartite_matching(A)
  prec <- TP / K_hat
  rec  <- TP / K_star
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, TP = TP, recall = rec, precision = prec, f1 = f1)
}

# Legacy bidirectional tagging metrics (for continuity checks)
tagging_metrics <- function(sel, causal, R, tau = 0.5) {
  K_hat <- length(sel)
  if (K_hat == 0L)
    return(data.frame(K_hat = 0L, recall = 0, precision = 0, f1 = 0))
  rec  <- mean(vapply(causal, function(cj) any(abs(R[sel, cj]) >= tau), logical(1)))
  prec <- mean(vapply(sel,   function(s)  any(abs(R[s, causal]) >= tau), logical(1)))
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, recall = rec, precision = prec, f1 = f1)
}
