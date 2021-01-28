#' Assigns taxonomy classifications to a provided table of sequences.
#'
#' @param seqtab The provided table of sequences. Can be generated from DADA2_main().
#' @param pathToDatabase A filepath to a genus-level classification database compatible with DADA2.
#' @param species A boolean specifying whether specie-level classification should be attempted.
#' @param pathToSpecies A filepath to a species-level classification database compatible with DADA2.
#'
#' @return Returns a table of classified sequences.
#' @export

label_taxonomy <- function(seqtab, pathToDatabase, species=FALSE, pathToSpecies=FALSE) {
    message("Assigning taxonomy...")
    taxa <- dada2::assignTaxonomy(seqtab.nochim, pathToDatabase, multithread=TRUE)

    if (species) {
        message("Assigning species...")
        taxa <- dada2::addSpecies(taxa, pathToSpecies)
    }

    return(taxa)
    message("Finished taxonomy annotation")
}
