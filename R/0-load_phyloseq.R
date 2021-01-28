#' A helper function to quickly load phyloseq objects. Currently unused.
#'
#' @param output_dir A filepath to the output directory, without the final slash.
#' @param rarefy_samples A boolean showing whether a rarefied phyloseq object should be loaded or not.
#'
#' @return The loaded phyloseq object.
#'
#' @keywords internal

loadphyloseq <- function(output_dir, rarefy_samples=FALSE) {
    if (rarefy_samples) {
        phyloseqobj <- readRDS(paste0(output_dir, fileplatform("/Rarefaction/RDS_objects/rarefied_ps.RDS")))
    } else {
        phyloseqobj <- readRDS(paste0(output_dir, fileplatform("/dada2/RDS_objects/phyloseq_output.RDS")))
    }
    return(phyloseqobj)
}
