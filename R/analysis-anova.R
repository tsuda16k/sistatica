#' function name
#'
#' function description
#'
#' @param data dataframe
#' @param col col name
#'
#' @return list
#' @export
#'
#' @examples
#' summarize_column(mtcars, "mpg")
summarize_column <- function(data, col) {
  if (!col %in% names(data)) {
    stop("error")
  }

  x <- data[[col]]

  list(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}
