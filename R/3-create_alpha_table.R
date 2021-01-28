#' Creates a table of alpha diversity metrics per sample. Optionally saves the table as a .tsv file.
#'
#' @param phyloseq_object The phyloseq object containing the target dataset.
#' @param outputpath (optional) A filepath to the folder that the table should be saved into.
#'
#' @return A table of alpha diversity metrics per sample.
#' @export

create_alpha_table <- function(phyloseq_object, outputpath=FALSE) {

    message("Starting microbiome diversity table generation...")

    tab <- microbiome::alpha(phyloseq_object, index = "all")

    if (outputpath) {
        dir.create(paste0(outputpath, fileplatform("/alpha_diversity_table")), showWarnings = FALSE)
        write.table(tab,
                    file=paste0(outputpath, fileplatform("/alpha_diversity_table/alpha_diversity.tsv")),
                    row.names=TRUE,
                    col.names=TRUE,
                    na = "",
                    sep = "\t",
                    quote = FALSE)
    }

    return(tab)
    message("Finished table generation\n")
}
