#' Generate Pastel Colors
#'
#' This function it's just for stetic purposes 😊
#'
#' @param n Number of colors we require
#' @return A list of n pastel colors.
#' @export
generate_pastel_colors <- function(n) {
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("The 'grDevices' package is required.")
  }

  # Generate pastel colors using HCL color space
  pastel_colors <- grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[-1], # Hue range
    c = 70,  # Chroma (controls color intensity)
    l = 85   # Lightness (higher for pastel effect)
  )

  return(pastel_colors)
}
