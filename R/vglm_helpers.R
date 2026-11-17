# Fit the vglm engine of mostm(). The model is built as in VGAM::vglm(), with
# two differences. A `constraints` list may name single columns of an inline
# registered rms basis term such as `rcs(time, 4)` or `lsp(time, 5)`
# (`"rcs(time, 4)time'"`), so each basis column can have its own constraint
# matrix in a partial proportional odds model; the fit then gives VGAM one
# assignment per basis column. Otherwise the fit uses VGAM's usual per-term
# assignment, because VGAM matches `parallel` formulas and `constraints` names
# against the assignment names, and a per-column assignment would silently
# ignore `parallel = FALSE ~ rcs(time, 4)` (#42). For the same reason a
# `parallel` formula naming a basis term whose single columns `constraints`
# names is refused rather than ignored. The returned fit always
# carries the per-column assignment, so prediction code sees the same object
# either way. The rms helpers `rcs()`, `lsp()`, and `%ia%` are available to
# the formula without `library(rms)`, so `time %ia% yprev` works as in rms.
#
# `setup` holds the arguments that mostm() has already resolved; the other
# formals are the VGAM::vglm() arguments a user may pass through mostm().
mostm_fit_vglm <- function(
  setup,
  family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
  na.action,
  etastart = NULL,
  mustart = NULL,
  coefstart = NULL,
  control = VGAM::vglm.control(...),
  model = FALSE,
  contrasts = NULL,
  constraints = NULL,
  extra = list(),
  form2 = NULL,
  qr.arg = TRUE,
  smart = TRUE,
  ...
) {
  data <- setup$data
  subset_rows <- setup$subset_rows
  function.name <- "vglm"

  formula <- add_rms_formula_helpers(setup$formula)
  if (!is.null(form2)) {
    form2 <- add_rms_formula_helpers(form2)
  }

  setup_smart <- utils::getFromNamespace("setup.smart", "VGAM")
  get_smart_prediction <- utils::getFromNamespace(
    "get.smart.prediction",
    "VGAM"
  )
  wrapup_smart <- utils::getFromNamespace("wrapup.smart", "VGAM")
  shadowvglm <- utils::getFromNamespace("shadowvglm", "VGAM")
  eval_vcontrol <- make_vgam_vcontrol_eval()
  get_xlevels <- utils::getFromNamespace(".getXlevels", "stats")

  smart_open <- FALSE
  if (smart) {
    setup_smart("write")
    smart_open <- TRUE
    on.exit(
      {
        if (smart_open) {
          wrapup_smart()
        }
      },
      add = TRUE
    )
  }

  # The row selection and weights were evaluated by mostm(). `na.action`,
  # `etastart`, and `mustart` are added to the model frame as in VGAM::vglm().
  mf <- list(quote(stats::model.frame), formula = formula, data = data)
  if (!is.null(subset_rows)) {
    mf$subset <- subset_rows
  }
  if (!is.null(setup$prior_weight_column)) {
    mf$weights <- as.name(setup$prior_weight_column)
  }
  for (nm in c("na.action", "etastart", "mustart")) {
    if (!is.null(setup$dots[[nm]])) {
      mf[[nm]] <- setup$dots[[nm]]
    }
  }
  mf$drop.unused.levels <- TRUE
  mf <- eval(as.call(mf), setup$eval_env)
  fit_data <- markov_align_model_data(data, mf)

  # Take the starting values from the model frame, as stats::glm() does, so
  # that they keep only the rows left after `subset` and `na.action`. Passing
  # the supplied values on unchanged would give them more rows than the fit.
  # Without `etastart` or `mustart` the frame has no such column and both stay
  # `NULL`.
  etastart <- stats::model.extract(mf, "etastart")
  mustart <- stats::model.extract(mf, "mustart")

  mt <- attr(mf, "terms")
  if (terms_has_offset(mt) || !is.null(stats::model.offset(mf))) {
    stop_unsupported_offset()
  }

  xlev <- get_xlevels(mt, mf)
  y <- stats::model.response(mf, "any")
  x <- if (!stats::is.empty.model(mt)) {
    stats::model.matrix(mt, mf, contrasts)
  } else {
    matrix(, NROW(y), 0)
  }

  # VGAM's cm.VGAM() and process.constraints() match the term labels of a
  # `parallel` formula and the names of a `constraints` list against
  # names(attr(x, "assign")). With one entry per basis column,
  # `parallel = FALSE ~ rcs(time, 4)` matches nothing and silently fits a
  # proportional odds spline, and a `constraints` list naming the whole term
  # fails with "'constraints' is half-pie" (#42). So the fit uses VGAM's usual
  # per-term assignment unless `constraints` names a single basis column, which
  # only the split assignment can express.
  split_assign <- split_rms_basis_assign(x, mt)
  constraints_name_basis_columns <- is.list(constraints) &&
    any(names(constraints) %in% split_assign$basis_columns)
  # With the split assignment a `parallel` formula naming a whole basis term
  # matches no column; mostm() refuses that combination once the family is
  # resolved below (see stop_parallel_formula_on_split_basis()).
  attr(x, "assign") <- if (constraints_name_basis_columns) {
    split_assign$assign
  } else {
    split_assign$term_assign
  }
  attr(x, "orig.assign.lm") <- split_assign$orig_assign

  if (!is.null(form2)) {
    if (!is.null(subset_rows)) {
      mostr_abort(
        "argument 'subset' cannot be used when argument 'form2' is used"
      )
    }
    retlist <- shadowvglm(
      formula = form2,
      family = family,
      data = data,
      na.action = na.action,
      control = VGAM::vglm.control(...),
      method = "vglm.fit",
      model = model,
      x.arg = TRUE,
      y.arg = TRUE,
      contrasts = contrasts,
      constraints = constraints,
      extra = extra,
      qr.arg = qr.arg
    )
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2
    if (length(Ym2) && NROW(Ym2) != NROW(y)) {
      mostr_abort("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2) && NROW(Xm2) != NROW(x)) {
      mostr_abort("number of rows of 'x' and 'Xm2' are unequal")
    }
  } else {
    Xm2 <- Ym2 <- NULL
  }

  offset <- 0

  w <- stats::model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  } else if (NCOL(w) == 1 && any(w < 0)) {
    mostr_abort("negative weights not allowed")
  }

  if (is.character(family)) {
    family <- get(family, envir = asNamespace("VGAM"), inherits = TRUE)
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "vglmff")) {
    mostr_abort("'family = ", family, "' is not a VGAM family function")
  }
  # VGAM's cm.VGAM() would match such a formula's basis term against the
  # per-column names, find nothing, and silently replace or keep the user's
  # column matrices against the formula's intent (#42).
  if (constraints_name_basis_columns) {
    stop_parallel_formula_on_split_basis(family, split_assign$basis_terms)
  }

  control <- eval_vcontrol(control, family, function.name, ...)
  if (length(methods::slot(family, "first"))) {
    eval(methods::slot(family, "first"))
  }

  fit <- VGAM::vglm.fit(
    x = x,
    y = y,
    w = w,
    offset = offset,
    Xm2 = Xm2,
    Ym2 = Ym2,
    etastart = etastart,
    mustart = mustart,
    coefstart = coefstart,
    family = family,
    control = control,
    constraints = constraints,
    extra = extra,
    qr.arg = qr.arg,
    Terms = mt,
    function.name = function.name,
    ...
  )

  # The returned fit carries the split assignment whichever one the fit used,
  # so `@assign`, `@x`, `@misc$orig.assign` and
  # VGAM::constraints(type = "term") list one entry per basis column. This is
  # safe because `fit$constraints` already holds one matrix per column of `x`,
  # and SOP prediction rebuilds the linear predictor from those matrices.
  attr(fit$x, "assign") <- split_assign$assign
  fit$misc$orig.assign <- split_assign$assign

  fit$misc$dataname <- setup$data_name
  if (smart) {
    fit$smart.prediction <- get_smart_prediction()
    wrapup_smart()
    smart_open <- FALSE
  }

  # mostm() replaces this call with its own.
  answer <- methods::new(
    Class = "vglm",
    assign = split_assign$assign,
    call = match.call(),
    coefficients = fit$coefficients,
    constraints = fit$constraints,
    criterion = fit$crit.list,
    df.residual = fit$df.residual,
    df.total = fit$df.total,
    dispersion = 1,
    effects = fit$effects,
    family = fit$family,
    misc = fit$misc,
    model = if (model) mf else data.frame(),
    R = fit$R,
    rank = fit$rank,
    residuals = as.matrix(fit$residuals),
    ResSS = fit$ResSS,
    smart.prediction = as.list(fit$smart.prediction),
    terms = list(terms = mt)
  )

  if (!smart) {
    answer@smart.prediction <- list(smart.arg = FALSE)
  }
  if (qr.arg) {
    class(fit$qr) <- "list"
    methods::slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts"))) {
    methods::slot(answer, "contrasts") <- attr(x, "contrasts")
  }
  if (length(fit$fitted.values)) {
    methods::slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  }
  methods::slot(answer, "na.action") <- if (
    length(aaa <- attr(mf, "na.action"))
  ) {
    list(aaa)
  } else {
    list()
  }
  methods::slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$weights)) {
    methods::slot(answer, "weights") <- as.matrix(fit$weights)
  }
  methods::slot(answer, "x") <- fit$x
  if (length(Xm2)) {
    methods::slot(answer, "Xm2") <- Xm2
  }
  if (length(Ym2)) {
    methods::slot(answer, "Ym2") <- as.matrix(Ym2)
  }
  if (!is.null(form2)) {
    methods::slot(answer, "callXm2") <- retlist$call
    answer@misc$Terms2 <- retlist$Terms2
  }
  answer@misc$formula <- formula
  answer@misc$form2 <- form2
  if (length(xlev)) {
    methods::slot(answer, "xlevels") <- xlev
  }
  methods::slot(answer, "y") <- as.matrix(fit$y)
  methods::slot(answer, "control") <- fit$control
  methods::slot(answer, "extra") <- if (length(fit$extra)) {
    if (is.list(fit$extra)) {
      fit$extra
    } else {
      warning("'extra' is not a list, therefore placing 'extra' into a list")
      list(fit$extra)
    }
  } else {
    list()
  }
  methods::slot(answer, "iter") <- fit$iter
  methods::slot(answer, "post") <- fit$post
  fit$predictors <- as.matrix(fit$predictors)
  if (length(fit$misc$predictors.names) == ncol(fit$predictors)) {
    dimnames(fit$predictors) <- list(
      dimnames(fit$predictors)[[1]],
      fit$misc$predictors.names
    )
  }
  methods::slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights)) {
    methods::slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)
  }

  attr(answer, "markov_vglm") <- TRUE
  attr(answer, "markov_split_assign") <- split_assign$has_basis

  list(fit = answer, fit_data = fit_data)
}

add_rms_formula_helpers <- function(formula) {
  if (!inherits(formula, "formula")) {
    return(formula)
  }

  formula_env <- environment(formula)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  helper_env <- new.env(parent = formula_env)

  for (handler in rms_basis_registry()) {
    helper <- handler$helper
    if (!exists(helper, envir = formula_env, inherits = TRUE)) {
      assign(
        helper,
        utils::getFromNamespace(helper, "rms"),
        envir = helper_env
      )
    }
  }
  if (!exists("%ia%", envir = formula_env, inherits = TRUE)) {
    assign(
      "%ia%",
      utils::getFromNamespace("%ia%", "rms"),
      envir = helper_env
    )
  }

  environment(formula) <- helper_env
  formula
}

# Returns two VGAM assignments of the columns of `x`: `term_assign`, one entry
# per model term as VGAM::vglm() builds it, and `assign`, where each inline rms
# basis term is replaced by one entry per generated column, named by the column
# (`"rms::rcs(time, 4)time'"`). `basis_columns` lists the names that only the
# split assignment has, which is how mostm_fit_vglm() recognises a
# `constraints` list that names single basis columns.
split_rms_basis_assign <- function(x, terms) {
  attrassigndefault <- utils::getFromNamespace("attrassigndefault", "VGAM")
  orig_assign <- attr(x, "assign")
  term_assign <- attrassigndefault(x, terms)
  assign <- term_assign

  basis_terms <- names(assign)[vapply(
    names(assign),
    function(term) length(rms_basis_handlers_for_term(term)) > 0L,
    logical(1)
  )]

  for (term in basis_terms) {
    cols <- assign[[term]]
    pos <- match(term, names(assign))
    assign[term] <- NULL
    split_cols <- stats::setNames(as.list(cols), colnames(x)[cols])
    assign <- append(assign, split_cols, after = pos - 1)
  }

  list(
    assign = assign,
    term_assign = term_assign,
    basis_terms = basis_terms,
    basis_columns = setdiff(names(assign), names(term_assign)),
    orig_assign = orig_assign,
    has_basis = length(basis_terms) > 0L
  )
}

# Term labels of the `parallel` formulas in a VGAM family. Families such as
# VGAM::cumulative() substitute their `parallel` argument into the
# `cm.VGAM(..., bool = <parallel>)` calls of the family's `constraints` slot,
# either as a formula object or as an unevaluated `~` call. The walk collects
# the labels of every formula it finds there; a logical `parallel` adds none.
vglm_family_parallel_term_labels <- function(family) {
  walk <- function(node) {
    if (inherits(node, "formula")) {
      return(attr(stats::terms(node), "term.labels"))
    }
    if (is.call(node) && identical(node[[1]], as.name("~"))) {
      return(attr(stats::terms(stats::as.formula(node)), "term.labels"))
    }
    if (is.call(node) || is.expression(node) || is.list(node)) {
      return(unlist(lapply(as.list(node), walk), use.names = FALSE))
    }
    character()
  }
  unique(walk(methods::slot(family, "constraints")))
}

stop_parallel_formula_on_split_basis <- function(family, basis_terms) {
  named <- intersect(vglm_family_parallel_term_labels(family), basis_terms)
  if (!length(named)) {
    return(invisible())
  }
  terms_text <- paste0("`", named, "`", collapse = ", ")
  mostr_abort(
    "`constraints` names single columns of the basis ",
    if (length(named) == 1L) "term " else "terms ",
    terms_text,
    ", and the family's `parallel` formula also names ",
    if (length(named) == 1L) "this term" else "these terms",
    ". A `parallel` formula treats a basis term as a whole, so it cannot be ",
    "combined with constraint matrices for its single columns. Give every ",
    "column of ",
    terms_text,
    " its constraint matrix in `constraints`, and use `parallel = FALSE` or a ",
    "`parallel` formula that does not name ",
    if (length(named) == 1L) "the term." else "these terms."
  )
}

make_vgam_vcontrol_eval <- function() {
  expr <- utils::getFromNamespace("vcontrol.expression", "VGAM")
  eval(
    substitute(
      function(control, family, function.name, ...) {
        eval(EXPR)
        control
      },
      list(EXPR = expr)
    ),
    envir = asNamespace("VGAM")
  )
}
