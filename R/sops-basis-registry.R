# The rms numeric basis helpers whose terms the vglm engine splits into one
# design column per basis piece.

rms_basis_registry <- function() {
  list(
    rcs = list(
      helper = "rcs",
      pattern = "(^|:)\\s*(rms::)?rcs\\("
    ),
    lsp = list(
      helper = "lsp",
      pattern = "(^|:)\\s*(rms::)?lsp\\("
    )
  )
}

rms_basis_handlers_for_term <- function(term) {
  handlers <- rms_basis_registry()
  names(handlers)[vapply(
    handlers,
    function(handler) grepl(handler$pattern, term, perl = TRUE),
    logical(1)
  )]
}
