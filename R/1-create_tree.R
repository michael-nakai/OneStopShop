#' Creates a phylogenetic tree using phangorn and DECEIPHER.
#'
#' @param seqtab A table of sequences. Can be generated via the [DADA2_main()] function.
#'
#' @return A phylogenetic tree object, used to make a full phyloseq object.
#'
#' @seealso [DADA2_main()] for sequence table creation.
#' @export

create_tree <- function(seqtab) {
    message("Creating phylogenetic tree...")

    seqs <- dada2::getSequences(seqtab)
    names(seqs) <- seqs # This propagates to the tip labels of the tree
    alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA)
    phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
    dm <- phangorn::dist.ml(phang.align)
    treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
    fit = phangorn::pml(treeNJ, data=phang.align)
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0)) # This takes a while

    return(fitGTR)
    message("Finished phylogenetic tree creation\n")
}
