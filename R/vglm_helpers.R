# Fit the vglm engine of mostm(). The model is built as in VGAM::vglm(), with
# two differences that the SOP prediction code relies on. Inline registered
# rms basis terms such as `rcs(time, 4)` or `lsp(time, 5)` get one VGAM
# assignment per generated column, so each basis column can have its own
# constraint matrix in a partial proportional odds model. The rms helpers
# `rcs()`, `lsp()`, and `%ia%` are available to the formula without
# `library(rms)`, so `time %ia% yprev` works as in rms.
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
  # `etastart`, and `mustart` are evaluated as VGAM::vglm() evaluates them.
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

  split_assign <- split_rms_basis_assign(x, mt)
  attr(x, "assign") <- split_assign$assign
  attr(x, "orig.assign.lm") <- split_assign$orig_assign

  if (!is.null(form2)) {
    if (!is.null(subset_rows)) {
      stop(
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
      stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2) && NROW(Xm2) != NROW(x)) {
      stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  } else {
    Xm2 <- Ym2 <- NULL
  }

  offset <- 0

  w <- stats::model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  } else if (NCOL(w) == 1 && any(w < 0)) {
    stop("negative weights not allowed")
  }

  if (is.character(family)) {
    family <- get(family, envir = asNamespace("VGAM"), inherits = TRUE)
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
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

  fit$misc$dataname <- setup$data_name
  if (smart) {
    fit$smart.prediction <- get_smart_prediction()
    wrapup_smart()
    smart_open <- FALSE
  }

  # mostm() replaces this call with its own.
  answer <- methods::new(
    Class = "vglm",
    assign = attr(x, "assign"),
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
  attr(answer, "markov_basis_terms") <- split_assign$basis_terms

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

split_rms_basis_assign <- function(x, terms) {
  attrassigndefault <- utils::getFromNamespace("attrassigndefault", "VGAM")
  orig_assign <- attr(x, "assign")
  assign <- attrassigndefault(x, terms)

  basis_terms <- names(assign)[vapply(
    names(assign),
    function(term) length(rms_basis_handlers_for_term(term)) > 0L,
    logical(1)
  )]
  basis_metadata <- vector("list", length(basis_terms))
  names(basis_metadata) <- basis_terms

  for (term in basis_terms) {
    cols <- assign[[term]]
    basis_metadata[[term]] <- list(
      handlers = rms_basis_handlers_for_term(term),
      columns = colnames(x)[cols],
      column_indices = as.integer(cols)
    )
    pos <- match(term, names(assign))
    assign[term] <- NULL
    split_cols <- stats::setNames(as.list(cols), colnames(x)[cols])
    assign <- append(assign, split_cols, after = pos - 1)
  }

  list(
    assign = assign,
    orig_assign = orig_assign,
    has_basis = length(basis_terms) > 0L,
    basis_terms = basis_metadata
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
