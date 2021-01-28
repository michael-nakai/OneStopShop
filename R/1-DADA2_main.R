#' Contains all DADA2 preprocessing, denoising, merging, and chimera removal processes.
#'
#' @param fastq_folderpath The filepath to the directory containing the FASTQ files to process.
#' @param outputpath The filepath to the directory where all filtered sequences and summary data should be output.
#' @param trim_left The number of base pairs to trim from the start of all sequences (min 0).
#' @param truncF The base number that forward reads should truncate at.
#' @param truncR The base number that reverse reads should truncate at.
#' @param forwardReadPattern A filename identifier that identifies FASTQ files with forward reads (ex: "_F").
#' @param reverseReadPattern A filename identifier that identifies FASTQ files with reverse reads (ex: "_R").
#' @param max_EE_allowed A vector containing the maximum error rates allowed (default c(2,2)).
#' @param max_N_allowed The maximum amount of unknown bases labelled N to allow before discarding a read (default 0).
#' @param min_length_allowed The minimum length a read can be. If the read length < min_length_allowed, the read is discarded (default 0).
#' @param min_quality_allowed The minimum quality score that a base can have before the entire read is discarded (default 0).
#' @param remove_phiX A boolean signifying whether reads assigned to phiX should be removed (default TRUE).
#' @param verbose_output A boolean signifying whether the main DADA2 step should be verbose or silent (default TRUE).
#' @param pooling A boolean signifying whether DADA2 should pool samples together prior to sample inference (default FALSE).
#'
#' @return Returns a table of sequences, all of which are non-chimeric and have been denoised and merged.
#'
#' @export

DADA2_main <- function(fastq_folderpath, outputpath, trim_left, truncF, truncR, forwardReadPattern, reverseReadPattern,
                       max_EE_allowed=c(2,2), max_N_allowed=0, max_length_allowed=10000, min_length_allowed=0,
                       min_quality_allowed=0, remove_phiX=TRUE, verbose_output=TRUE, pooling=FALSE) {

    message("Starting file sort...")

    # Sort into forward and reverse reads
    fnFs <- sort(list.files(fastq_folderpath,
                            pattern=forwardReadPattern,
                            full.names = TRUE))

    fnRs <- sort(list.files(fastq_folderpath,
                            pattern=reverseReadPattern,
                            full.names = TRUE))

    # Extract sample names, assuming the sample names are before the first underscore
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

    # Make the output paths per filtered file, then label them with the corresponding sample names
    dir.create(paste0(outputpath, fileplatform("/dada2/")), showWarnings = FALSE)
    filtFs <- file.path(outputpath, fileplatform("dada2/filtered/"), paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(outputpath, fileplatform("dada2/filtered/"), paste0(sample.names, "_R_filt.fastq.gz"))

    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

    # Do the filtering and save the output
    message("Starting DADA2 filtering and trimming...")

    trunclens <- c(truncF, truncR)
    out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                         truncLen = trunclens,
                         maxN = max_N_allowed,
                         maxEE = max_EE_allowed,
                         maxLen = max_length_allowed,
                         minLen = min_length_allowed,
                         minQ = min_quality_allowed,
                         truncQ = 2,
                         rm.phix = remove_phiX,
                         compress = TRUE, multithread = TRUE, verbose = verbose_output)

    dada2_save_path <- paste0(outputpath, fileplatform("/dada2/summary/filterAndTrimOutput.csv"))
    dir.create(paste0(outputpath, fileplatform("/dada2/summary/")), showWarnings = FALSE)
    write.csv(out, dada2_save_path)

    # Check the estimated error rates for possible base --> base errors
    # This takes a while, depending on how many reads we have, so fire and forget it
    # If the red/black lines trend in the same direction and fit relatively well together, then everything's OK
    message("Checking estimated base --> base error rates...")
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    dada2::plotErrors(errF, nominalQ=TRUE)

    # Apply the core sample inference algorithm (again, can take a while)
    message("Applying sample inference algorithm...")
    dadaFs <- dada2::dada(filtFs, err=errF, multithread=TRUE, pool = pooling)
    dadaRs <- dada2::dada(filtRs, err=errR, multithread=TRUE, pool = pooling)

    # Merge the paired-end reads
    message("Merging forward and reverse reads...")
    mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=verbose_output)

    # Make a sequence table
    message("Making sequence table...")
    seqtab <- dada2::makeSequenceTable(mergers)

    # Check that many of the sequences are all around the same length
    table_save_path <- paste0(outputpath, fileplatform("/dada2/summary/readLength.csv"))
    lengthTable <- table(nchar(getSequences(seqtab)))
    write.csv(lengthTable, table_save_path)

    # Remove chimeras
    message("Removing chimeras...")
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=verbose_output)

    # Output summary of all steps, then save
    message("Creating DADA2 summary...")
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    dada2_final_save_path <- paste0(outputpath, fileplatform("/dada2/summary/Summary.csv"))
    write.csv(track, dada2_final_save_path)

    # Save the final seqtab
    dir.create(paste0(outputpath, fileplatform("/dada2/RDS_objects/")), showWarnings = FALSE)
    saveRDS(seqtab.nochim, file = paste0(outputpath, fileplatform("/dada2/RDS_objects/seqtab_nochim.RDS")))

    return(seqtab.nochim)
    message("Finished DADA2 denoising\n")
}
