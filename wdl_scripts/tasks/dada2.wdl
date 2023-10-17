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
        File dada2_classifier
    }

    command<<<
        set -ex -o pipefail

        R --no-save <<Rscript

        library(dada2)
        ##https://benjjneb.github.io/dada2/tutorial.html

        fnFs <- "~{fastq_1}"
        fnRs <- "~{fastq_2}"
        sample.names <- "~{sample_name}"

        # Place filtered files in filtered/ subdirectory
        filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
        filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))

        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        trimLeft=c(~{trim_f},~{trim_r}),
                        truncLen=c(~{trunc_f},~{trunc_r}),
                        maxN=0, maxEE=c(2,5), truncQ=~{trunc_q}, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) #отрезать концы с низким качеством прочтения

        errF <- learnErrors(filtFs,multithread = TRUE)
        errR <- learnErrors(filtRs,multithread = TRUE)

        derepFs <- derepFastq(filtFs, verbose=FALSE)
        derepRs <- derepFastq(filtRs, verbose=FALSE)

        # Name the derep-class objects by the sample names
        names(derepFs) <- sample.names
        names(derepRs) <- sample.names

        dadaFs <- dada(derepFs, err=errF, multithread=TRUE,verbose = FALSE,pool="pseudo") #https://benjjneb.github.io/dada2/pseudo.html про псевдопулинг
        dadaRs <- dada(derepRs, err=errR, multithread=TRUE,verbose = FALSE,pool="pseudo") #псевдопулинг увеличивает врем¤ (примерно в 2-2,5 раза)

        mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

        seqtab <- makeSequenceTable(mergers)

        seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)

        taxa <- assignTaxonomy(seqtab.nochim, "~{dada2_classifier}", multithread=TRUE, verbose=FALSE)

        ##add species by exact match
        #taxa2 <- addSpecies(taxa, "~/Documents/R/dada2/silva_species_assignment_v138.1.fa.gz")

        taxa.print <- taxa # Removing sequence rownames for display only
        rownames(taxa.print) <- NULL

        ##remove seqs from rownames
        seqtab.print<-t(seqtab.nochim)
        rownames(seqtab.print) <- NULL

        #output
        write.table(data.frame(seqtab.print, taxa.print), file='seqtab.nochim.tsv', quote=FALSE, sep='\t', col.names = NA)

        ##stat
        getN <- function(x) sum(getUniques(x))
        if(length(sample.names) >1) {
        track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        } else {
          track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
        }
        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track) <- sample.names
        write.table(track, file='stat.tsv', quote=FALSE, sep='\t', col.names = NA)

        ##delete filtered temp dir
        unlink("filtered", recursive = TRUE)

        Rscript
    >>>

  runtime {
      docker: "blekhmanlab/dada2:latest"
  }

  output {
     File stat_tsv = "stat.tsv"
     File seqtab_nochim_tsv = "seqtab.nochim.tsv"
  }
}


