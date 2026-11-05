# Allowed-transition matrix, transport feasibility, and IPF fit diagnostics.
#
# The allowed-transition matrix Z (spatial rulebook, Principle 6, item 4) marks class
# pairs that can never receive flow. It is passed as a table of pairs; pairs not
# listed are allowed, staying in a class is always allowed, and exempt classes
# (balancing classes such as NODATA) keep all their transitions. Every argument is
# optional: without a matrix, nothing changes.

#' Normalize an allowed-transition table
#'
#' @param allowed `data.frame` with columns `from_class`, `to_class` and `allowed`
#'   (0/1 or logical). Pairs not listed are allowed.
#' @param exempt_classes Classes that keep all their transitions (e.g. balancing classes).
#' @return A `data.table` of the disallowed pairs (`from_class`, `to_class`), excluding
#'   stays and exempt classes; `NULL` if `allowed` is `NULL`.
#' @export
normalize_allowed_transitions <- function(allowed, exempt_classes = "NODATA") {
  if (is.null(allowed)) return(NULL)
  allowed <- data.table::as.data.table(allowed)
  need <- c("from_class", "to_class", "allowed")
  if (!all(need %in% names(allowed))) {
    stop("allowed transitions need columns ", paste(need, collapse = ", "), ".")
  }
  a <- allowed$allowed
  if (is.logical(a)) a <- as.integer(a)
  if (any(is.na(a)) || !all(a %in% c(0, 1))) stop("allowed transitions: 'allowed' must be 0/1 or TRUE/FALSE.")
  dup <- allowed[, .N, by = .(from_class, to_class)][N > 1]
  if (nrow(dup)) stop("allowed transitions: duplicated pairs, e.g. ", dup$from_class[1], " -> ", dup$to_class[1], ".")
  dis <- unique(allowed[a == 0, .(from_class = as.character(from_class), to_class = as.character(to_class))])
  stays <- dis[from_class == to_class]
  if (nrow(stays)) message("allowed transitions: staying in a class is always allowed; ignoring ", nrow(stays), " stay row(s).")
  dis[from_class != to_class & !(from_class %in% exempt_classes) & !(to_class %in% exempt_classes)]
}

# TRUE for each (from, to) pair that is disallowed.
.pair_disallowed <- function(from, to, disallowed) {
  if (is.null(disallowed) || !nrow(disallowed)) return(rep(FALSE, length(from)))
  paste(from, to, sep = "\r") %in% paste(disallowed$from_class, disallowed$to_class, sep = "\r")
}

#' Transport feasibility of two margins on a support
#'
#' Checks whether flows on the cells where `support` is positive can meet the row
#' totals `supply` and the column totals `demand` (a bipartite max-flow problem).
#'
#' @param support Matrix (rows = from classes, columns = to classes); positive cells may carry flow.
#' @param supply,demand Non-negative row and column totals with equal sums.
#' @param tol Relative tolerance on the shortfall.
#' @return A list with `feasible`, `max_flow`, `total`, and `shortfall`.
#' @export
check_transport_feasibility <- function(support, supply, demand, tol = 1e-9) {
  supply <- as.numeric(supply); demand <- as.numeric(demand)
  n <- length(supply); k <- length(demand)
  if (!all(dim(support) == c(n, k))) stop("support must be length(supply) x length(demand).")
  total <- sum(supply)
  if (total <= 0) return(list(feasible = TRUE, max_flow = 0, total = 0, shortfall = 0))
  # Nodes: 1 = source, 2..(n+1) = rows, (n+2)..(n+k+1) = columns, n+k+2 = sink
  N <- n + k + 2; src <- 1; snk <- N
  cap <- matrix(0, N, N)
  cap[src, 1 + seq_len(n)] <- supply
  cap[1 + n + seq_len(k), snk] <- demand
  big <- 2 * total
  cap[1 + seq_len(n), 1 + n + seq_len(k)] <- ifelse(support > 0, big, 0)
  eps <- tol * total
  flow <- 0
  repeat {
    parent <- rep(0L, N); parent[src] <- src
    queue <- src
    while (length(queue) && parent[snk] == 0L) {
      u <- queue[1]; queue <- queue[-1]
      nxt <- which(cap[u, ] > eps & parent == 0L)
      parent[nxt] <- u
      queue <- c(queue, nxt)
    }
    if (parent[snk] == 0L) break
    path_cap <- Inf; v <- snk
    while (v != src) { u <- parent[v]; path_cap <- min(path_cap, cap[u, v]); v <- u }
    v <- snk
    while (v != src) { u <- parent[v]; cap[u, v] <- cap[u, v] - path_cap; cap[v, u] <- cap[v, u] + path_cap; v <- u }
    flow <- flow + path_cap
  }
  shortfall <- total - flow
  list(feasible = shortfall <= eps, max_flow = flow, total = total, shortfall = shortfall)
}

# Session registry of solve_gross_transitions() fits.
.ipf_registry <- new.env(parent = emptyenv())
.ipf_registry$log <- NULL
.ipf_registry$warned <- 0L

.record_ipf_fit <- function(row) {
  .ipf_registry$log <- data.table::rbindlist(list(.ipf_registry$log, row), use.names = TRUE, fill = TRUE)
  invisible(NULL)
}

#' Diagnostics of the gross-transition fits in this session
#'
#' One row per `solve_gross_transitions()` call: number of classes, feasibility,
#' iterations, stop reason, remaining relative margin error, and whether the fit met
#' its margins (`converged`). Fits that miss their margins are also warned about
#' (the first 10 per session).
#'
#' @param reset If `TRUE`, clear the registry (and the warning counter) after returning it.
#' @return A `data.table` (empty if no fit was recorded).
#' @export
ipf_diagnostics <- function(reset = FALSE) {
  out <- if (is.null(.ipf_registry$log)) data.table::data.table() else data.table::copy(.ipf_registry$log)
  if (isTRUE(reset)) { .ipf_registry$log <- NULL; .ipf_registry$warned <- 0L }
  out
}
