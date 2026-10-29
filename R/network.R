#' Define a Node of a Network
#'
#' A node is one sub-process of a DMU. It consumes inputs and/or intermediate
#' products and produces final outputs and/or intermediate products. Variables
#' are referred to by the column names of the data frame passed to [ndea()].
#' An intermediate product links two nodes by appearing as `z_out` in the node
#' that produces it and as `z_in` in the node that consumes it.
#'
#' @param name Node name (a short label).
#' @param inputs Column names of the final inputs the node uses.
#' @param outputs Column names of the final outputs the node produces.
#' @param z_in Column names of the intermediate products the node consumes.
#' @param z_out Column names of the intermediate products the node produces.
#'
#' @return An object of class `"ndea_node"`.
#' @seealso [network()], [ndea()]
#' @examples
#' node("A", inputs = c("labour", "capital"), z_out = "interim")
#' node("B", z_in = "interim", outputs = "sales")
#' @export
node <- function(name, inputs = character(), outputs = character(),
                 z_in = character(), z_out = character()) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a single non-empty string.", call. = FALSE)
  }
  nd <- list(name = name, inputs = as.character(inputs),
             outputs = as.character(outputs), z_in = as.character(z_in),
             z_out = as.character(z_out))
  if (length(nd$inputs) + length(nd$z_in) == 0L) {
    stop("Node '", name, "' needs at least one input or incoming intermediate (`inputs` or `z_in`).",
         call. = FALSE)
  }
  if (length(nd$outputs) + length(nd$z_out) == 0L) {
    stop("Node '", name, "' needs at least one output or outgoing intermediate (`outputs` or `z_out`).",
         call. = FALSE)
  }
  structure(nd, class = "ndea_node")
}

#' @export
print.ndea_node <- function(x, ...) {
  f <- function(v) if (length(v)) paste(v, collapse = ", ") else "-"
  cat("Node '", x$name, "': inputs [", f(x$inputs), "], z_in [", f(x$z_in),
      "] -> outputs [", f(x$outputs), "], z_out [", f(x$z_out), "]\n", sep = "")
  invisible(x)
}

#' Define a Network of Nodes
#'
#' Collects [node()] definitions into a network structure for [ndea()].
#'
#' Checks performed: node names are unique; a final input or output belongs to
#' only one node (shared inputs are not supported); an intermediate that is only
#' produced or only consumed gives a warning.
#'
#' @param ... [node()] objects.
#'
#' @return An object of class `"ndea_network"`.
#' @seealso [node()], [ndea()]
#' @examples
#' net <- network(
#'   node("A", inputs = c("labour", "capital"), z_out = "interim"),
#'   node("B", z_in = "interim", outputs = "sales")
#' )
#' net
#' @export
network <- function(...) {
  nodes <- list(...)
  if (length(nodes) < 2L) stop("A network needs at least two nodes.", call. = FALSE)
  if (!all(vapply(nodes, inherits, logical(1), "ndea_node"))) {
    stop("All arguments of `network()` must be created with `node()`.", call. = FALSE)
  }
  nms <- vapply(nodes, function(nd) nd$name, character(1))
  if (anyDuplicated(nms)) stop("Node names must be unique.", call. = FALSE)
  names(nodes) <- nms
  fin <- unlist(lapply(nodes, function(nd) c(nd$inputs, nd$outputs)))
  if (anyDuplicated(fin)) {
    stop("Final input/output variable(s) used by more than one node: ",
         paste(unique(fin[duplicated(fin)]), collapse = ", "),
         ". Shared inputs are not supported.", call. = FALSE)
  }
  zin <- unlist(lapply(nodes, function(nd) nd$z_in))
  zout <- unlist(lapply(nodes, function(nd) nd$z_out))
  if (anyDuplicated(zout)) {
    stop("An intermediate is produced by more than one node: ",
         paste(unique(zout[duplicated(zout)]), collapse = ", "), call. = FALSE)
  }
  both_final <- intersect(fin, c(zin, zout))
  if (length(both_final)) {
    stop("Variable(s) used both as final and as intermediate: ",
         paste(both_final, collapse = ", "), call. = FALSE)
  }
  lonely <- c(setdiff(zout, zin), setdiff(zin, zout))
  if (length(lonely)) {
    warning("Intermediate(s) not both produced and consumed within the network: ",
            paste(lonely, collapse = ", "), call. = FALSE)
  }
  structure(list(nodes = nodes), class = "ndea_network")
}

#' @export
print.ndea_network <- function(x, ...) {
  cat("Network with", length(x$nodes), "nodes\n")
  for (nd in x$nodes) print(nd)
  invisible(x)
}
