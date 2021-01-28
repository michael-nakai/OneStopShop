#' Plots all general alpha diversity metrics on one plot.
#'
#' @param phyloseq_object The phyloseq object containing the data to be visualized.
#'
#' @return A plot showing all general metrics of alpha diversity.
#' @export

plot_alpha_all <- function(phyloseq_object) {
    general_richness <- phyloseq::plot_richness(phyloseq_object)
    return(general_richness)
}
