knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.8,
  out.width = "100%",
  message = FALSE,
  warning = FALSE
)

# Output longer than 10 lines is collapsed, unless the chunk sets
# `fold_output = FALSE`.
local({
  output_hook <- knitr::knit_hooks$get("output")
  knitr::knit_hooks$set(output = function(x, options) {
    output <- output_hook(x, options)
    if (
      !knitr::is_html_output() ||
        isFALSE(options$fold_output) ||
        length(strsplit(x, "\n")[[1]]) <= 10L
    ) {
      return(output)
    }
    paste0(
      "\n\n<details><summary>Show output</summary>\n\n",
      output,
      "\n\n</details>\n\n"
    )
  })
})

# Chunks with `fold = TRUE` show their code collapsed; output stays visible.
local({
  source_hook <- knitr::knit_hooks$get("source")
  knitr::knit_hooks$set(source = function(x, options) {
    source <- source_hook(x, options)
    if (!isTRUE(options$fold) || !knitr::is_html_output()) {
      return(source)
    }
    paste0(
      "\n\n<details><summary>Show code</summary>\n\n",
      source,
      "\n\n</details>\n\n"
    )
  })
})
