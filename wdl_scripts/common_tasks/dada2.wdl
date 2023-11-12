version 1.0

task Dada2 {

    input {
        File? fastq_1
        File? fastq_2
        String sample_name
        Int trunc_q
        Int trim_f
        Int trim_r
        Int trunc_f
        Int trunc_r
        Int minOverlap
        File dada2_classifier
        String docker
    }

    command<<<
        set -ex -o pipefail

        R --no-save <<Rscript

        library(dada2)
        ##https://benjjneb.github.io/dada2/tutorial.html

        fnFs <- "~{fastq_1}"
        fnRs <- "~{fastq_2}"

        # Place filtered files in filtered/ subdirectory
        filtFs <- file.path("filtered", paste0("~{sample_name}", "_F_filt.fastq.gz"))
        filtRs <- file.path("filtered", paste0("~{sample_name}", "_R_filt.fastq.gz"))

        # Filter and trim
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        trimLeft=c(~{trim_f},~{trim_r}),
                        truncLen=c(~{trunc_f},~{trunc_r}),
                        maxN=0, maxEE=c(2,5), truncQ=~{trunc_q}, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)

        # If the filtered list is NOT empty
        if (any(out[, "reads.out"] != 0)) {

          # Learn the Error Rates
          errF <- learnErrors(filtFs,multithread = TRUE)
          errR <- learnErrors(filtRs,multithread = TRUE)

          # In new versions dereplication is performed "on the fly"
          # But it still be a useful diagnostic to look at how many unique sequences were in total reads.
          #derepFs <- derepFastq(filtFs, verbose=FALSE)
          #derepRs <- derepFastq(filtRs, verbose=FALSE)

          dadaFs <- dada(filtFs, err=errF, multithread=TRUE,verbose = FALSE,pool="pseudo") #https://benjjneb.github.io/dada2/pseudo.html про псевдопулинг
          dadaRs <- dada(filtRs, err=errR, multithread=TRUE,verbose = FALSE,pool="pseudo") #псевдопулинг увеличивает врем¤ (примерно в 2-2,5 раза)

          # Merge paired reads
          mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, maxMismatch=5,
                               minOverlap='~(minOverlap}', verbose=TRUE)

          # If a majority of reads failed to merge, you may need to revisit the truncLen parameter
          # And make sure that the truncated reads span your amplicon. Or reduce minOverlap.
          # Non-overlapping reads are supported, but not recommended with justConcatenate=TRUE
          if(all(mergers$accept == FALSE)) {
              mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                                    justConcatenate=TRUE, verbose=TRUE)
          }

          # Construct sequence table
          seqtab <- makeSequenceTable(mergers)

          # Remove chimeras
          seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

          # Assign taxonomy
          taxa <- assignTaxonomy(seqtab.nochim, "~{dada2_classifier}", multithread=TRUE, verbose=FALSE)

          ##add species by exact match
          #taxa2 <- addSpecies(taxa, "~/Documents/R/dada2/silva_species_assignment_v138.1.fa.gz")

          taxa.print <- taxa # Removing sequence rownames for display only
          rownames(taxa.print) <- NULL

          # remove seqs from rownames
          seqtab.print<-t(seqtab.nochim)
          rownames(seqtab.print) <- NULL

          # output
          write.table(data.frame(seqtab.print, taxa.print), file='~{sample_name}_seqtab.nochim.tsv', quote=FALSE, sep='\t', col.names = NA)

          # stat
          getN <- function(x) sum(getUniques(x))
          track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
          colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
          rownames(track) <- "~{sample_name}"
          write.table(track, file='~{sample_name}_stat.tsv', quote=FALSE, sep='\t', col.names = NA)

        } else {
          message <- "In filterAndTrim no reads passed the filter. Please revisit your filtering parameters."
          write.table(message, '~{sample_name}_stat.tsv', sep = "\t", col.names = FALSE, quote = FALSE)
        }

        # delete filtered temp dir
        unlink("filtered", recursive = TRUE)

        Rscript
    >>>

  runtime {
      docker: "blekhmanlab/dada2:latest"
  }

  output {
     File stat_tsv = "~{sample_name}_stat.tsv"
     File? seqtab_nochim_tsv = "~{sample_name}_seqtab.nochim.tsv"
  }
}
