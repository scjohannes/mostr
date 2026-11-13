# User-facing documentation of the package options. The defaults live where
# each option is read: `getOption("mostr.<name>", <default>)`. Keep the values
# below in step with those calls.

#' Options that limit memory use
#'
#' mostr splits large calculations into pieces so that memory stays bounded.
#' These options set the size of the pieces and the limits at which mostr
#' changes method or stops. The defaults suit most analyses; change them only
#' when a calculation is too slow or runs out of memory. The first two options
#' change only how the work is split, never the results. The others stop a
#' calculation, or move it to a slower method, when it would need more memory
#' than the limit.
#'
#' Set an option with [options()] before the call and restore it afterwards,
#' as in the example. Sizes in bytes are easiest to write as multiples of
#' `1024^2` (one MiB).
#'
#' @section Options that change only how the work is split:
#' \describe{
#'   \item{`mostr.max_batched_design_rows`}{Default `5e4` rows. For each visit,
#'     the state occupancy probabilities need one design row for every profile
#'     and every state a patient can be in before the visit. Up to this many
#'     rows per visit, mostr builds them in one batch; above it, one previous
#'     state at a time, which is faster and needs less memory for large
#'     problems. The results are identical either way.}
#'   \item{`mostr.second_order_working_bytes`}{Default `32 * 1024^2` (32 MiB).
#'     The working memory per piece in the R engine's second-order recursion,
#'     which works through the pairs of states at the last two visits a few at
#'     a time. Peak memory is about six times this value. Larger values give
#'     no measurable speed-up. The results are identical for any value.}
#' }
#'
#' @section Options that stop a calculation or change its method:
#' \describe{
#'   \item{`mostr.execution_plan_max_bytes`}{Default `256 * 1024^2` (256 MiB).
#'     The compiled engine keeps the design matrices of all visits in memory.
#'     When they would need more than this, predictions use the R engine
#'     instead, which gives the same results more slowly (about 3.6 times
#'     slower for MVN inference in one measurement), and delta-method
#'     inference stops with an error.}
#'   \item{`mostr.delta_max_bytes`}{Default `256 * 1024^2` (256 MiB). The
#'     largest single block of memory that delta-method inference
#'     ([inferences()] with `method = "delta"`) may request: the derivatives
#'     of each estimate with respect to each model coefficient, their
#'     workspaces, and covariance blocks requested with [stats::vcov()]. A
#'     calculation that would need more stops with an error that gives the
#'     size it needs. Requesting fewer rows, for example fewer `times`, or
#'     averaging with [avg_sops()] instead of keeping one row per profile,
#'     reduces this memory.}
#'   \item{`mostr.posterior_working_bytes`}{Default `256 * 1024^2` (256 MiB).
#'     For `rmsb::blrm()` fits, posterior draws are processed in groups whose
#'     working memory stays below this value. Smaller values use less memory
#'     and, in measurements, cost no time.}
#'   \item{`mostr.max_sops_draw_cells`}{Default `5e7` values, about 400 MB.
#'     With `rmsb::blrm()` fits, [sops()] and [comparisons()] keep one value
#'     per posterior draw, profile, visit and state. A call that would keep
#'     more stops with an error. Use fewer `n_draws`, or average with
#'     [avg_sops()] or [avg_comparisons()], before raising it.}
#' }
#'
#' These limits count the largest numeric arrays only. The fitted model, the
#' prediction data, temporary copies made by R and other overhead need memory
#' as well, so a calculation can use more memory than any one limit even when
#' every check passes. Raise a limit only after checking that the reported
#' size fits comfortably in the available memory.
#'
#' @name mostr_options
#' @aliases mostr.max_batched_design_rows
#' @aliases mostr.second_order_working_bytes
#' @aliases mostr.execution_plan_max_bytes
#' @aliases mostr.delta_max_bytes
#' @aliases mostr.posterior_working_bytes
#' @aliases mostr.max_sops_draw_cells
#' @seealso [inferences()], [sops()], [avg_sops()]
#' @examples
#' # Allow delta-method inference to use up to 512 MiB, then restore the
#' # previous setting.
#' old_options <- options(mostr.delta_max_bytes = 512 * 1024^2)
#' options(old_options)
NULL
