knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.8,
  out.width = "100%",
  message = FALSE,
  warning = FALSE
)

local({
  output_hook <- knitr::knit_hooks$get("output")
  knitr::knit_hooks$set(output = function(x, options) {
    output <- output_hook(x, options)
    if (!knitr::is_html_output() || length(strsplit(x, "\n")[[1]]) <= 10L) {
      return(output)
    }
    paste0(
      "\n\n<details><summary>Show output</summary>\n\n",
      output,
      "\n\n</details>\n\n"
    )
  })
})
