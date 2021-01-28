#' Rarefies samples in the provided phyloseq object to the specified rarefaction depth. Optionally saves the resulting phyloseq object.
#'
#' @param phyloseq_object The phyloseq object to rarefy.
#' @param rarefaction_depth A number specifying the number of reads to rarefy to per sample.
#' @param outputpath (optional) A filepath specifying the folder to save the rarefied object to.
#'
#' @return A rarefied phyloseq object.
#' @export

rarefy_phyloseq <- function(phyloseq_object, rarefaction_depth, outputpath=FALSE) {

    message("Starting rarefaction...")

    # Rarefy the phyloseq object and save it
    ps.rarefied = rarefy_even_depth(phyloseq_object, sample.size=rarefaction_depth, replace=FALSE)

    if (outputpath) {
        # Create the output folder
        dir.create(paste0(outputpath, fileplatform("/Rarefaction")), showWarnings = FALSE)
        dir.create(paste0(outputpath, fileplatform("/Rarefaction/RDS_objects")), showWarnings = FALSE)
        saveRDS(ps.rarefied, file = paste0(outputpath, fileplatform("/Rarefaction/RDS_objects/rarefied_ps.RDS")))
    }

    return(ps.rarefied)
    message("Finished rarefying samples")
}
