#' Creates a list of alpha diversity plots based on the metrics specified.
#'
#' @param phyloseq_object The phyloseq object containing the relevant dataset.
#' @param alpha_types A vector containing the types of alpha diversity metrics to show.
#'
#' @return Returns a list of alpha diversity plots, or raises an error if alpha_types isn't a vector
#'
#' @seealso <https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness> for the primarily used function.
#' @export

plot_alpha_specific <- function(phyloseq_object, alpha_types=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) {

    # Check that alpha_types is a vector
    if (!(is.vector(alpha_types))) {
        stop("alpha_types is not a vector!")
    }

    outputList <- vector("list", length = length(alpha_types))
    i <- 1

    for (divType in alpha_types) {
        outputList[[i]] <- phyloseq::plot_richness(ps_object, measures = divType)
        i <- i + 1
    }

    return(outputList)
}
