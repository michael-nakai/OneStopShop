#' Creates a phyloseq object from the supplied data. If outputpath is specified, save the phyloseq object.
#'
#' @param sequence_table A sequence table to include. Can be generated from [DADA2_main()].
#' @param metadata A metadata table to include. Can be generated from [format_q2_metadata()].
#' @param taxa_table A taxonomy table to include. Can be generated from [label_taxonomy()].
#' @param tree (optional) A phylogenetic tree to include. Can be generated from [create_tree()].
#' @param outputpath (optional) A filepath to a folder to output the saved phyloseq object as an RDS object.
#' @param substitute_ASV_names (optional) A boolean specifying whether ASV names should be substituted with ASV1, ASV2... The names will be stored in "dna".
#'
#' @return A phyloseq object consisting of the supplied sequence table, metadata, taxonomy table, and (optional) phylogenetic tree.
#'
#' @seealso [DADA2_main()] for sequence table creation, [format_q2_metadata()] for metadata importing, [label_taxonomy()] for taxonomy table creation, and [create_tree()] for phylogenetic tree creation.
#' @export

handoff_to_phyloseq <- function(sequence_table, metadata, taxa_table, tree=FALSE, outputpath=FALSE, substitute_ASV_names=FALSE){

    message("Starting phyloseq handoff...")

    if (tree == FALSE) {
        ps_object <- phyloseq::phyloseq(phyloseq::otu_table(seqtab, taxa_are_rows=FALSE),
                                        phyloseq::sample_data(metadata),
                                        phyloseq::tax_table(taxa_table))
    } else {
        ps_object <- phyloseq::phyloseq(phyloseq::otu_table(seqtab, taxa_are_rows=FALSE),
                                        phyloseq::sample_data(metadata),
                                        phyloseq::phy_tree(tree),
                                        phyloseq::tax_table(taxa_table))
    }

    if (substitute_ASV_names) {
        dna <- Biostrings::DNAStringSet(taxa_names(ps_object))
        names(dna) <- phyloseq::taxa_names(ps_object)
        ps_object <- phyloseq::merge_phyloseq(ps_object, dna)
        phyloseq::taxa_names(ps_object) <- paste0("ASV", seq(phyloseq::ntaxa(ps_object)))
    }

    if (outputpath) {
        dir.create(paste0(outputpath, fileplatform("/dada2/")), showWarnings = FALSE)
        dir.create(paste0(outputpath, fileplatform("/dada2/RDS_objects/")), showWarnings = FALSE)
        saveRDS(ps_object, file = paste0(outputpath, fileplatform("/dada2/RDS_objects/phyloseq_output.RDS")))
    }

    return(ps_object)
    message("Finished phyloseq handoff\n")
}
