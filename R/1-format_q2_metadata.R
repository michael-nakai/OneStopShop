#' Formats metadata from a QIIME2-compatible .tsv file to a data frame, and prepares it for downstream use.
#'
#' @param metadata_path A filepath to the metadata file.
#'
#' @return A dataframe of the imported metadata.
#' @export

format_q2_metadata <- function(metadata_path) {
    # Format the metadata for the phyloseq sample_data
    message("Formatting metadata for phyloseq handoff...")

    temp_sample_data <- as.data.frame(fread(metadata_path))
    if (temp_sample_data[1,1] == "#q2:types") {   # If the second row is #q2:types then delete it
        temp_sample_data <- temp_sample_data[-c(1),]
    }
    rowname_sample_data <- temp_sample_data[,-1]  # Reassign the first column as the rownames
    rownames(rowname_sample_data) <- temp_sample_data[,1]

    return(rowname_sample_data)
    message("Finished formatting metadata\n")
}
