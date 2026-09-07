# Run from the repository root with Rscript .github/check-vignette-output.R.
local({
  path <- tempfile(fileext = ".Rmd")
  writeLines(
    c(
      "---",
      'title: "Output folding check"',
      "output: rmarkdown::html_vignette",
      "vignette: >",
      "  %\\VignetteIndexEntry{Output folding check}",
      "---",
      "```{r, include=FALSE}",
      'source("_setup.R", local = TRUE)',
      "```",
      "```{r}",
      'cat("Short result\\n")',
      'cat(paste0("Long result ", 1:11, "\\n"), sep = "")',
      "plot(1:3)",
      "```"
    ),
    path
  )
  html <- rmarkdown::render(
    path,
    knit_root_dir = normalizePath("vignettes"),
    quiet = TRUE
  )
  page <- xml2::read_html(html)
  details <- xml2::xml_find_all(page, "//details")
  stopifnot(
    length(details) == 1L,
    is.na(xml2::xml_attr(details, "open")),
    trimws(xml2::xml_text(xml2::xml_find_first(details, "summary"))) ==
      "Show output",
    grepl("Long result 11", xml2::xml_text(details), fixed = TRUE),
    !grepl("Short result", xml2::xml_text(details), fixed = TRUE),
    length(xml2::xml_find_all(page, "//details//img")) == 0L,
    length(xml2::xml_find_all(page, "//img")) == 1L,
    length(xml2::xml_find_all(page, "//details//pre[contains(@class, 'r')]")) ==
      0L
  )
  cat(
    "Long output is expandable; short output, code, and plots stay visible.\n"
  )
})
