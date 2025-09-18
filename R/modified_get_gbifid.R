# Load the taxize library first
library(taxize)

#' A patched version of taxize::get_gbifid that correctly passes '...'
#'
#' The standard taxize::get_gbifid function does not pass the ... argument to
#' its internal helper, preventing low-level options from being used. This
#' patched version corrects that single oversight, allowing the pass-through
#' of arguments like `opts` to fix HTTP/2 connection errors.
#'
#' @param sci (character) One or more scientific names.
#' @param ask (logical) Run in interactive mode.
#' @param ... Additional arguments passed down to the internal crul client.
#'   Use `opts = list(http_version = 2)` here.
#' @param sciname Deprecated, see `sci`.
#' @param messages (logical) Print messages.
#' @param rows (numeric) Number of rows to consider.
#' @param phylum,class,order,family,rank Optional taxonomic filters.
#' @param method (character) one of "backbone" or "lookup".
#'
#' @return A list of 'gbifid' class objects.
get_gbifid_patched <- function(sci, ask = TRUE, ..., sciname = NULL, messages = TRUE,
                               rows = NA, phylum = NULL, class = NULL, order = NULL,
                               family = NULL, rank = NULL, method = "backbone") {
  
  # --- This code is copied directly from taxize::get_gbifid ---
  
  checkmate::assert(sci, "character")
  if (!is.null(sciname)) {
    lifecycle::deprecate_warn(
      when = "v0.9.97",
      what = "get_gbifid(sciname)",
      with = "get_gbifid(sci)"
    )
    sci <- sciname
  }
  
  # This internal function `fxn` is where the original bug is.
  # We are fixing it by adding `...` to the `taxize:::get_gbifid_` call.
  fxn <- function(x, ...) {
    taxize:::get_gbifid_( # Note the ::: to access the internal function
      x,
      ask = ask,
      messages = messages,
      rows = rows,
      phylum = phylum,
      class = class,
      order = order,
      family = family,
      rank = rank,
      method = method,
      # --- THIS IS THE ONLY CHANGE ---
      # The original function was missing this line, which is why opts was ignored.
      ...
    )
  }
  
  # --- The rest of the original function code is unchanged ---
  att <- attributes(sci)
  if (inherits(sci, "taxon_state")) {
    out <- fxn(sci, ...)
  } else {
    out <- lapply(sci, fxn, ...) # lapply passes the ... from the main call to fxn
    names(out) <- sci
    class(out) <- "gbifid"
  }
  
  if (!is.null(att)) {
    att <- att[names(att) %in% c("match", "multiple_matches", "pattern_match")]
    attributes(out) <- c(attributes(out), att)
  }
  
  return(out)
}
