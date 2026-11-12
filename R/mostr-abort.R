# Raising errors.
#
# Every error that mostr raises goes through mostr_abort(), never stop(). The
# rule is in AGENTS.md; "Errors" in ARCHITECTURE.md describes the design.

# Raise an error of class `mostr_error`.
#
# `...` are message parts, pasted together without a separator as stop() does,
# so `mostr_abort("`x` must have length ", n, ".")` reads like the stop() call
# it replaced. The message is used literally: braces and `%` are not
# interpolated, and there are no cli bullets.
#
# `class` adds classes in front of `mostr_error`, so callers can catch one kind
# of error. `data` is a named list of extra fields stored in the condition, read
# back with `condition$field`.
#
# `call` is the call shown as "Error in `f()`:". The default,
# mostr_error_call(), is the call of the exported function the user called, so
# users never see the name of an internal function. `NULL` shows no call.
#
# A single error condition, as in `mostr_abort(e)` inside an error handler, is
# re-raised unchanged, gaining class `mostr_error` if it lacks it.
mostr_abort <- function(
  ...,
  class = NULL,
  data = list(),
  call = mostr_error_call()
) {
  if (...length() == 1L && inherits(..1, "error")) {
    condition <- ..1
    if (!inherits(condition, "mostr_error")) {
      base_class <- intersect(class(condition), c("error", "condition"))
      class(condition) <- c(
        setdiff(class(condition), base_class),
        "mostr_error",
        base_class
      )
    }
    rlang::cnd_signal(condition)
  }
  force(call)
  # Fields are spliced into rlang::abort(), so a field named like one of its
  # arguments (`body`, `parent`, ...) would silently change the error instead.
  fields <- names(data) %||% rep("", length(data))
  if (any(fields == "") || any(fields %in% names(formals(rlang::abort)))) {
    rlang::abort(
      "`data` must hold named fields that aren't rlang::abort() arguments.",
      .internal = TRUE
    )
  }
  message <- paste(unlist(lapply(list(...), as.character)), collapse = "")
  rlang::abort(
    message,
    class = c(class, "mostr_error"),
    !!!data,
    call = call,
    use_cli_format = FALSE
  )
}

# The call of the exported mostr function the user called, or NULL.
#
# Walks the call stack for frames that run an exported function. It takes the
# innermost one called from outside the package: that is the call the user
# wrote, even when `sops(fit)` in `inferences(sops(fit))` is evaluated inside
# another exported function. If every such frame was called from package code
# (test files run in an environment inside the namespace, for example), it
# takes the outermost one. With no exported function on the stack, as in a
# parallel worker or a direct call of an internal function, it returns NULL and
# no call is shown.
#
# A call made through do.call() holds the function itself rather than its name,
# which rlang would not show, so the exported name replaces it.
mostr_error_call <- function() {
  namespace <- topenv(environment(mostr_error_call))
  exported <- mget(
    getNamespaceExports(namespace),
    envir = namespace,
    mode = "function",
    ifnotfound = list(NULL),
    inherits = FALSE
  )
  exported <- Filter(Negate(is.null), exported)
  in_namespace <- function(fn) {
    is.function(fn) &&
      !is.primitive(fn) &&
      identical(topenv(environment(fn)), namespace)
  }
  export_index <- function(fn) {
    if (!in_namespace(fn)) {
      return(NA_integer_)
    }
    match(TRUE, vapply(exported, identical, logical(1), fn))
  }

  # The last frame is this function's own.
  n_frames <- sys.nframe() - 1L
  parents <- sys.parents()
  outermost <- NULL
  chosen <- NULL
  for (i in seq_len(n_frames)) {
    index <- export_index(sys.function(i))
    if (is.na(index)) {
      next
    }
    if (is.null(outermost)) {
      outermost <- c(i, index)
    }
    parent <- parents[[i]]
    if (parent == 0L || !in_namespace(sys.function(parent))) {
      chosen <- c(i, index)
    }
  }
  found <- chosen %||% outermost
  if (is.null(found)) {
    return(NULL)
  }
  call <- sys.call(found[[1L]])
  head <- call[[1L]]
  named <- is.symbol(head) ||
    (is.call(head) && identical(head[[1L]], as.name("::")))
  if (!named) {
    call[[1L]] <- as.name(names(exported)[[found[[2L]]]])
  }
  call
}
