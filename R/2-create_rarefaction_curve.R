#' Creates a rarefaction curve using vegan's rarecurve().
#'
#' @param phyloseq_object The phyloseq object to be loaded.
#' @param step_size The intervals the rarefaction curve should be crafted on.
#'
#' @return A rarefaction curve.
#'
#' @seealso <https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/rarefy>
#' @export

create_rarefaction_curve <- function(phyloseq_object, step_size=20) {

    message("Starting rarefaction curve generation...")

    rarefaction_curve <- vegan::rarecurve(otu_table(ps_object, FALSE), step=step_size, cex=0.5)

    return(rarefaction_curve)
    message("Finished rarefaction curve generation\n")
}
