#' Filter biological data by variance, or find a cutoff
#'
#' @description This function calculates a cutoff based on data variance, aiming to filter out the "0 variance peak".
#' @param data a numeric matrix or a numeric vector. If matrix - Genes must be in rows
#' @param mads numeric (0-10]. Number of MADs to be added to the median for cutoff. Default=2 according to the outlier detection threshold. 
#' @param must_genes a string vector with gene names that must be included in the output. Relevant only for matrix input.
#' @param plot_density logical. Whether to plot the density of data variance distribution with the calculated cutoffs.
#' @param ... additional arguments to other functions
#'
#' @return If `plot_density` == FALSE: a filtered matrix (if `data` is a matrix) or a numeric cutoff (if `data` is a vector); If `plot_density` == TRUE: a `ggplot` object.
#'
#' @export
#' @importFrom modeest asselin
#' @importFrom stats median density mad var
#' @importFrom assertthat assert_that
#' @import ggplot2
#'
#' @examples
#'\dontrun{
#' madvar(df, mads = 1, FALSE, must_genes = c("ERBB2", "BRCA1"), plot_density = TRUE)
#'}
madvar <- function (data, mads = 2, must_genes = NULL, plot_density = FALSE, ...) {
  assert_that(mads > 0 && mads <= 10, msg = "Error: 'mads' must be > 0 and <= 10")
  abline <- axis <- par <- text <- NULL
  mat <- FALSE
  if (is.matrix(data) | is.data.frame(data))
    mat <- TRUE
  if (mat) {
    gene_variance <- apply(as.matrix(data), 1, var, na.rm = T)
  }
  else {
    gene_variance <- data
  }
  med <- median(gene_variance)
  if (med <= 0) {
    warning("Most features have 0 variance.\nI'm removing them and re-calculating the cutoff.")
    if (mat) {
      data <- data[gene_variance > 0, ]
      gene_variance <- apply(as.matrix(data), 1, var, na.rm = T)
    }
    else {
      gene_variance <- gene_variance[gene_variance > 0]
    }
  }
  MAD <- mad(gene_variance)
  median_peak <- 2 * med - min(gene_variance)
  if (plot_density) { # Exploratoration mode with plot
    assert_that(mat, msg = "Error: The plot can be generated only for a matrix/data frame")
    ggplot(as.data.frame(gene_variance), aes(x=gene_variance)) +
      geom_density(alpha=0.4) +
      labs(x = "Gene variance", y = "Density", ...) +
      theme_bw() +
      geom_vline(xintercept = med, color = "gray", linetype = "dashed") +
      geom_vline(xintercept = asselin(gene_variance, ...), color = "blue", linetype = "dashed") +
      geom_vline(xintercept = med + MAD * mads, color = "red", linetype = "dashed")
  } else {
    if (mat) {
      data[gene_variance >= med + MAD * mads | rownames(data) %in% must_genes, ]
    } else {
      med + MAD * mads
    }
  }
}
