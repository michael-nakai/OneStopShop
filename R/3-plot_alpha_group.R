#' Creates a list of alpha diversity plots with samples separated into groups based on the provided grouping variable.
#'
#' @param phyloseq_object The phyloseq object containing the dataset to be plotted.
#' @param grouping_variable The column name in the metadata that the samples should be grouped by. Multiple columns can be looped over if provided in a vector.
#' @param alpha_types The types of alpha diversity metrics to plot.
#'
#' @return Returns a list of alpha diversity plots (if grouping variable isn't a vector) or a list of lists (if grouping variable is a vector). Raises an error if alpha_types isn't a vector
#' @export

plot_alpha_group <- function(phyloseq_object, grouping_variable, alpha_types=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) {

    # Check that alpha_types is a vector
    if (!(is.vector(alpha_types))) {
        stop("alpha_types is not a vector!")
    }

    if (!(is.vector(grouping_variable))) {

        outputList <- vector("list", length = 0)

        for (divType in alpha_types) {
            outputList[[divType]] <- phyloseq::plot_richness(ps_object, x = grouping_variable, measures = divType)
        }
    }

    if (is.vector(grouping_variable)) {
        outputList <- vector("list", length = 0)

        for (variab in grouping_variable) {
            intermList <- vector("list", length = 0)

            for (divType in alpha_types) {
                intermList[[divType]] <- phyloseq::plot_richness(ps_object, x = variab, measures = divType)
            }

            outputList[[variab]] <- intermList
        }
    }

    return(outputList)
}
