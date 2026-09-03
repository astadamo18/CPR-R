# Shared lm()-like formula/data resolution for cpr() and pcpr().
#
# Deliberately uses all.vars() rather than model.frame()/model.matrix(): the
# integrated regressors, stationary regressors, and deterministic terms here
# are plain numeric series, not factors needing contrasts or interactions,
# so full formula semantics would be more machinery than the problem needs.
# One consequence worth stating plainly: a transformed term like log(x1) in
# a formula is *not* evaluated -- all.vars(log(x1) ~ ...) just returns "x1",
# so a formula RHS must name existing columns of `data` directly. Create the
# transformed column in `data` first (e.g. `df$log_x1 <- log(df$x1)`) and
# reference that name instead.

#' Extract `y` and `x` from a two-sided formula and a data frame
#'
#' @param formula A two-sided formula, `y ~ x1 + x2`. The left-hand side
#'   must be a single column; the right-hand side names one or more columns
#'   of `data` verbatim (see file-level comment on transformed terms).
#' @param data A data frame containing the named columns.
#' @return A list with `y` (vector) and `x` (matrix, columns named after
#'   the right-hand-side variables, in the order written).
#' @keywords internal
extract_formula_xy <- function(formula, data) {
  if (is.null(data)) {
    stop("`data` must be supplied when the model is given as a formula ",
         "(e.g. cpr(y ~ x1 + x2, data = df, orders = 2)).", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (length(formula) != 3) {
    stop("The model formula must be two-sided, e.g. `y ~ x1 + x2`.", call. = FALSE)
  }

  resp_name <- all.vars(formula[[2]])
  if (length(resp_name) != 1) {
    stop("The left-hand side of the formula must be a single column, ",
         "e.g. `y ~ x1 + x2`.", call. = FALSE)
  }
  rhs_vars <- all.vars(formula[[3]])
  if (length(rhs_vars) == 0) {
    stop("The right-hand side of the formula must name at least one column, ",
         "e.g. `y ~ x1 + x2`.", call. = FALSE)
  }

  missing_cols <- setdiff(c(resp_name, rhs_vars), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  y <- data[[resp_name]]
  x <- as.matrix(data[rhs_vars])
  colnames(x) <- rhs_vars

  list(y = y, x = x)
}

#' Resolve a one-sided formula's variables against a data frame into a matrix
#'
#' Used for `w`/`deter` when given as e.g. `~ w1 + w2` instead of a raw
#' matrix.
#' @keywords internal
resolve_formula_vars <- function(formula, data, label) {
  if (is.null(data)) {
    stop("`data` must be supplied when `", label, "` is given as a formula.", call. = FALSE)
  }
  vars <- all.vars(formula)
  if (length(vars) == 0) {
    stop("The `", label, "` formula must name at least one column, e.g. `~ ", label, "1`.",
         call. = FALSE)
  }
  missing_cols <- setdiff(vars, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) named in `", label, "` not found in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  out <- as.matrix(data[vars])
  colnames(out) <- vars
  out
}

#' Resolve a single column reference (`id`/`time` for pcpr()) against `data`
#'
#' If `ref` is a length-1 character string naming a column of `data`, return
#' that column; otherwise return `ref` unchanged (already a raw vector, the
#' original calling convention).
#' @keywords internal
resolve_column_ref <- function(ref, data) {
  if (!is.null(data) && !is.null(ref) && is.character(ref) && length(ref) == 1 && ref %in% names(data)) {
    return(data[[ref]])
  }
  ref
}

#' Resolve a multi-column reference (`w`/`deter` for pcpr()) against `data`
#'
#' If `ref` is a character vector whose entries all name columns of `data`,
#' return them as a matrix; otherwise return `ref` unchanged (already a raw
#' matrix/vector, the original calling convention).
#' @keywords internal
resolve_columns_ref <- function(ref, data) {
  if (!is.null(data) && !is.null(ref) && is.character(ref) && all(ref %in% names(data))) {
    out <- as.matrix(data[ref])
    colnames(out) <- ref
    return(out)
  }
  ref
}
