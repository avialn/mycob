version 1.0

task FastQC {

    input {
        File? fastq
        String docker
    }

    command <<<
        set -euxo pipefail

        chmod 755 /opt/fastqc/fastqc
        JAVA=$(which java)

        /opt/fastqc/fastqc \
        ~{fastq} \
        -java=$JAVA \
        --noextract \
        --outdir "."
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File qc_report_html = select_first(glob("*.html"))
    }
}

task Trimmomatic {

    input {
        File fastq_1
        File fastq_2
        Int threads
        String docker
    }

    command <<<
        set -euxo pipefail

        java -jar /opt/trimmomatic/trimmomatic-0.39.jar PE \
        -threads ~{threads} \
        ~{fastq_1} ~{fastq_2} \
        R1_trimmed.fastq.gz R1_unpaired.fastq.gz \
        R2_trimmed.fastq.gz R2_unpaired.fastq.gz \
        ILLUMINACLIP:/opt/trimmomatic/adapters.fa:2:30:5 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File trim_fastq_1 = "R1_trimmed.fastq.gz"
        File trim_fastq_2 = "R2_trimmed.fastq.gz"
    }
}

task Cutadapt {

    input {
        File fastq_1
        File fastq_2
        File? adapter_f
        File? adapter_r
        String docker
    }

    command <<<
        set -euxo pipefail

        . /opt/cutadapt/bin/activate

        /opt/cutadapt/bin/cutadapt \
        --cores=0 \
        --report=minimal \
        -a file:~{adapter_f} \
        -A file:~{adapter_r} \
        -o cut_R1.fastq.gz -p cut_R2.fastq.gz \
        ~{fastq_1} \
        ~{fastq_2} > cutadapt_summary.tsv
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File cut_fastq_1 = "cut_R1.fastq.gz"
        File cut_fastq_2 = "cut_R2.fastq.gz"
        File cutadapt_summary_tsv = "cutadapt_summary.tsv"
    }
}

task Bowtie2Filter {

    input {
        File? fastq_1
        File? fastq_2
        File index_tar
        Int threads
        String bowtie2_options = "--very-sensitive-local"
        String docker
    }

    command <<<
        set -euxo pipefail

        count="$(cat ~{fastq_1} ~{fastq_2} | wc -l)"
        count=$((count / 4))

        jq --null-input --arg count "$count" '{"row_reads":$count}' > 'row_reads.count'

        tar xf '~{index_tar}' -C /tmp

        /opt/bowtie2/bowtie2 \
        -x /tmp/GRCh38_ERCC/GRCh38_ERCC \
        --very-sensitive-local \
        -p ~{threads} \
        -1 ~{fastq_1} \
        -2 ~{fastq_2} \
        -S bowtie2.sam
    >>>

    output {
        File row_reads_count = "row_reads.count"
        File output_sam = "bowtie2.sam"
    }

    runtime {
        docker: "~{docker}"
    }
}

task SamtoolsFilter {

    input {
        File bowttie2_output_sam
        Int threads
        String docker
    }

    command <<<
        set -euxo pipefail

        samtools sort -n \
        -o "bowtie2_host.bam" \
        -@ ~{threads} -T /tmp ~{bowttie2_output_sam}

        # Extract reads [pairs] that did NOT map to the index
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 \
        -1 'host_filtered1.fastq' \
        -2 'host_filtered2.fastq' \
        -0 /dev/null -s /dev/null bowtie2_host.bam

        count="$(cat host_filtered{1,2}.fastq | wc -l)"
        count=$((count / 4))

        jq --null-input --arg count "$count" '{"host_filtered_reads":$count}' > 'host_filtered_reads.count'
    >>>

    output {
       File host_filtered_fastq_1 = "host_filtered1.fastq"
       File host_filtered_fastq_2 = "host_filtered2.fastq"
       File host_filtered_reads_count = "host_filtered_reads.count"
       File bowtie2_host_bam = "bowtie2_host.bam"
    }

    runtime {
      docker: "~{docker}"
    }
}

task RgiBwt {
    input {
        File fastq_1
        File fastq_2
        File card_json
        String docker
    }

    command <<<
        set -exuo pipefail

        rgi bwt -1 ~{fastq_1} \
        -2 ~{fastq_2} \
        -a kma \
        -o amr_report --clean
    >>>

    output {
        File kma_amr_results_txt = "amr_report.allele_mapping_data.txt"
        File artifacts_mapping_stats = "amr_report.artifacts_mapping_stats.txt"
        File gene_mapping_data = "amr_report.gene_mapping_data.txt"
        File overall_mapping_stats = "amr_report.overall_mapping_stats.txt"
        File reference_mapping_stats = "amr_report.reference_mapping_stats.txt"
        File output_sorted_length_100 = "amr_report.sorted.length_100.bam"
        File output_sorted_length_100_bai = "amr_report.sorted.length_100.bam.bai"
        File kma_amr_results_json = "amr_report.allele_mapping_data.json"
    }

    runtime {
        docker:"~{docker}"
    }
}

task Flash {
    input {
        File fastq_1
        File fastq_2
        String docker
    }

    command <<<
        set -exuo pipefail

        /opt/FLASH-1.2.11/flash \
        ~{fastq_1} ~{fastq_2} \
        -o flash_output
    >>>

    output {
        File output_fq = "flash_output.extendedFrags.fastq"
    }

    runtime {
        docker: "~{docker}"
    }
}

task SeqtkToFa {
    input {
        File flash_output_fq
        String docker
    }

    command <<<
        set -exuo pipefail

        /opt/seqtk/seqtk seq \
        -a ~{flash_output_fq} \
        > output.fa
    >>>

    output {
        File output_fa = "output.fa"
    }

    runtime {
        docker: "~{docker}"
    }
}

task Spades {

    input {
        File seqtk_output_fa
        String docker
    }

    command <<<
        set -euxo pipefail
        function handle_failure()
        {
            echo ";ASSEMBLY FAILED" > spades/contigs.fasta
            echo ";ASSEMBLY FAILED" > spades/scaffolds.fasta
            exit 0
        }
        trap handle_failure ERR

        python3 /opt/SPAdes-3.15.5-Linux/bin/spades.py \
        -s ~{seqtk_output_fa} \
        -o "spades/" \
        -m 100 -t 36 \
        --only-assembler 1>&2
    >>>

    output {
        File contigs_fa = "spades/contigs.fasta"
        File? scaffolds_fa = "spades/scaffolds.fasta"
    }

    runtime {
        docker: "~{docker}"
    }
}

task SeqtkFilter {
    input {
        File contigs
        Int min_contig_length
        String docker
    }

    command <<<
        set -euxo pipefail
        if [[ $(head -n 1 "~{contigs}") == ";ASSEMBLY FAILED" ]]; then
            echo ";ASSEMBLY FAILED" > contigs_filtered.fasta
        else
            /opt/seqtk/seqtk seq -L ~{min_contig_length} ~{contigs} > contigs_filtered.fasta
        fi
    >>>

    output {
        File contigs_filtered = "contigs_filtered.fasta"
    }

    runtime {
        docker: "~{docker}"
    }
}

task RgiMain {

    input {
        File contigs_filt
        File card_json
        String docker
    }

    command <<<
        set -exuo pipefail

        if [[ $(head -n 1 "~{contigs_filt}") == ";ASSEMBLY FAILED" ]]; then
            # simulate empty outputs
            echo "{}" > contig_amr_report.json
            cp /tmp/empty-main-header.txt contig_amr_report.txt
        else
            rgi main \
            -i "~{contigs_filt}" \
            -o contig_amr_report \
            -t contig -a BLAST --include_nudge --clean
        fi
    >>>

    output {
        Array[File] output_main = glob("contig_amr_report*")
        File output_json = "contig_amr_report.json"
        File output_txt = "contig_amr_report.txt"
    }

    runtime {
        docker: "~{docker}"
    }
}

task GeneCoverage {
    input {
        File main_amr_results
        File main_output_json
        String docker
    }

    command <<<
    set -euxo pipefail
    python3 <<CODE
    import pandas as pd
    import numpy as np
    import json
    df = pd.read_csv("~{main_amr_results}", delimiter="\t").loc[:, ["ORF_ID", "ID", "Model_ID", "Hit_Start", "Hit_End", "Percentage Length of Reference Sequence"]]

    # create seq length reference map
    with open("~{main_output_json}") as json_file:
        rgi_main_json = json.load(json_file)
    db_seq_length = {}
    for ind, row in df.iterrows():
        db_seq_length[row["Model_ID"]] = len(rgi_main_json[row["ORF_ID"]][row["ID"]]["dna_sequence_from_broadstreet"])

    agg_res = df.groupby(["ID", "Model_ID"]).agg(lambda x: list(x))

    gene_coverage = []
    for ind, row, in agg_res.iterrows():
      max_end = -1
      gene_coverage_bps = 0
      for start, end, in sorted(zip(row["Hit_Start"], row["Hit_End"]), key=lambda x: x[0]):
        gene_coverage_bps += max(0, # if end-max(max_end, start) is negative
                      # segment is contained in a previous segment, so don't add
                     end - max(max_end, start) # if max_end > start, don't double count the already added portion of the gene
                    )
        max_end = max(max_end, end)
      gene_coverage.append({
        "ID": ind[0],
        "gene_coverage_bps": gene_coverage_bps,
        "db_seq_length": db_seq_length[ind[1]],
        "gene_coverage_perc": np.round((gene_coverage_bps/db_seq_length[ind[1]])*100, 4)
      })
    if gene_coverage:
        gene_coverage_df = pd.DataFrame(gene_coverage)
        gene_coverage_df.to_csv("gene_coverage.tsv", index=None, sep="\t")
    else:
        gene_coverage_df = pd.DataFrame(columns=["ID", "gene_coverage_bps", "db_seq_length", "gene_coverage_perc"])
        gene_coverage_df.to_csv("gene_coverage.tsv", index=None, sep="\t")
    CODE
    >>>
    output {
        File output_tsv = "gene_coverage.tsv"
    }

    runtime {
        docker: "~{docker}"
    }
}

task Report {

    input {
        File main_output
        File kma_output
        File gene_coverage_tsv
        String sample_name
        String docker
    }

    command <<<
        set -euxo pipefail

        python3 <<CODE
        import pandas as pd
        pd.set_option('mode.chained_assignment', None)

        def clean_aro(df, column_name):
            """modifies dataframe inplace to clean the ARO string"""
            df[column_name] = df[column_name].map(lambda x: x.strip())


        def append_suffix_to_colname(df, suffix):
            """append suffix to column name in place"""
            df.rename(lambda x: x + suffix, axis="columns", inplace=True)


        def this_or_that(df, this_colname, that_colname):
            """returns this if not nan, else returns that"""
            thisorthat = df.apply(
                lambda x: x[this_colname] if not pd.isna(x[this_colname]) else x[that_colname],
                axis=1,
            )
            if thisorthat.empty:
                return ""
            return thisorthat


        main_output = pd.read_csv("~{main_output}", sep="\t")
        clean_aro(main_output, "Best_Hit_ARO")
        main_output.sort_values(by='Best_Hit_ARO', inplace=True)

        kma_output = pd.read_csv("~{kma_output}", sep="\t")
        clean_aro(kma_output, "ARO Term")
        kma_output.sort_values(by = 'ARO Term', inplace=True)

        # rename main-amr and main-species columns to account for duplciate column names
        append_suffix_to_colname(main_output, "_contig_amr")

        # vv get rid of weird additional space in some of the contig amr fields to make matching work downstream
        main_output["Contig_contig_amr"] = main_output["Contig_contig_amr"].map(
            lambda x: x.strip()
        )
        main_output["ARO_contig"] = main_output["Best_Hit_ARO_contig_amr"]

        main_output.to_csv("main_output.tsv", index=None, sep="\t")

        main_output['ARO_contig_amr'] = main_output['ARO_contig_amr'].astype(str)

        append_suffix_to_colname(kma_output, "_kma_amr")  # ARO Term

        kma_output['ARO Accession_kma_amr'] = kma_output['ARO Accession_kma_amr'].astype(str)

        # remove kma results from variant models (because these results are inaccurate)
        kma_output = kma_output[kma_output['Reference Model Type_kma_amr'] != "protein variant model"] # remove protein variant model
        kma_output = kma_output[kma_output['Reference Model Type_kma_amr'] != "rRNA gene variant model"] # remove rRNA variant model

        kma_output["ARO_kma"] = kma_output["ARO Term_kma_amr"]

        kma_output.to_csv("kma_output.tsv", index=None, sep="\t")

        # final merge of MAIN and KMA combined results
        merge = main_output.merge(kma_output, left_on = 'ARO_contig_amr',
                        right_on = 'ARO Accession_kma_amr', how='outer',
                        suffixes = [None, None])

        merge["ARO_accession"] = this_or_that(merge, "ARO_contig_amr", "ARO Accession_kma_amr")
        merge["ARO_overall"] = this_or_that(merge, "ARO_contig", "ARO_kma")
        merge["Gene_Family_overall"] = this_or_that(
        merge, "AMR Gene Family_contig_amr", "AMR Gene Family_kma_amr"
        )
        merge["Drug_Class_overall"] = this_or_that(
        merge, "Drug Class_contig_amr", "Drug Class_kma_amr"
        )
        merge["model_overall"] = this_or_that(
        merge, "Model_type_contig_amr", "Reference Model Type_kma_amr"
        )
        merge["Resistance_Mechanism_overall"] = this_or_that(
        merge, "Resistance Mechanism_contig_amr", "Resistance Mechanism_kma_amr"
        )
        merge.to_csv("comprehensive_AMR_metrics.tsv", index=None, sep='\t')

        df = pd.read_csv("comprehensive_AMR_metrics.tsv", sep='\t')

        big_table = df[['ARO_overall', 'Gene_Family_overall', 'Drug_Class_overall', 'Resistance_Mechanism_overall',
                'model_overall', 'All Mapped Reads_kma_amr', 'Percent Coverage_kma_amr','Depth_kma_amr',
                'Cut_Off_contig_amr', 'Percentage Length of Reference Sequence_contig_amr',
                'Best_Identities_contig_amr']]

        big_table.sort_values(by='Gene_Family_overall', inplace=True)

        big_table.to_csv("bigtable_report.tsv", sep='\t', index=None)

        if big_table.empty: # if no outputs, simply return
            open("primary_AMR_report.tsv", "a").close()
            exit(0)

        gene_coverage = pd.read_csv("~{gene_coverage_tsv}", sep='\t')
        def remove_na(input_set):
            set_list = list(input_set)
            return([i for i in set_list if i == i])

        this_list = list(set(df['ARO_overall']))
        result_df = {}
        for index in this_list:#[0:1]:
            sub_df = df[df['ARO_overall'] == index]
            result = {}
            gf = remove_na(set(sub_df['AMR Gene Family_contig_amr']).union(set(sub_df['AMR Gene Family_kma_amr'])))
            result['sample_name'] = "~{sample_name}"
            result['gene_family'] = ';'.join(gf) if len(gf) > 0 else None

            dc = remove_na(set(sub_df['Drug Class_contig_amr']).union(set(sub_df['Drug Class_kma_amr'])))
            result['drug_class'] = ';'.join(dc) if len(dc) > 0 else None
            rm = remove_na(set(sub_df['Resistance Mechanism_contig_amr']).union(set(sub_df['Resistance Mechanism_kma_amr'])))
            result['resistance_mechanism'] = ';'.join(rm) if len(rm) > 0 else None

            sub_df.loc[(sub_df['Cut_Off_contig_amr'] == 'Strict') & (sub_df['Nudged_contig_amr'] == True), 'Cut_Off_contig_amr'] = "Nudged"
            co = remove_na(set(sub_df['Cut_Off_contig_amr']))

            result['cutoff'] = ';'.join(co) if len(co) > 0 else None

            mt = remove_na(set(sub_df['Model_type_contig_amr']).union(set(sub_df['Reference Model Type_kma_amr'])))
            result['model_type'] = ';'.join([' '.join(i.split(' ')[0:-1]) for i in mt]) if len(mt) > 0 else None

            rcb = remove_na(set(sub_df['Percent Coverage_kma_amr']))
            result['read_coverage_breadth'] = max(rcb) if len(rcb) > 0 else None

            gene_id = ";".join(remove_na(set(sub_df["ID_contig_amr"].unique())))
            if gene_id:
                contig_coverage = gene_coverage[gene_coverage["ID"] == gene_id]["gene_coverage_perc"].iloc[0]
            else:
                contig_coverage = None
            result["contig_coverage_breadth"] = contig_coverage

            cd = remove_na(set(sub_df['Depth_kma_amr']))
            result['read_coverage_depth'] = max(cd) if len(cd) > 0 else None

            result['num_contigs'] = len(remove_na(set(sub_df['Contig_contig_amr'])))

            nr = remove_na(set(sub_df['All Mapped Reads_kma_amr']))
            result['num_reads'] = max(nr) if len(nr) > 0 else None

            read_gi = remove_na(set(sub_df["Reference Sequence_kma_amr"]))
            result["read_gene_id"] = ";".join(read_gi) if len(read_gi) > 0 else None

            pid = remove_na(set(sub_df['Best_Identities_contig_amr']))
            result['contig_percent_id'] = sum(pid) / len(pid) if len(pid) > 0 else None

            result_df[index] = result

        final_df = pd.DataFrame.from_dict(result_df)
        final_df = final_df.transpose()
        final_df = final_df[["sample_name", "gene_family", "drug_class", "resistance_mechanism",
                   "model_type", "num_contigs", "cutoff", "contig_coverage_breadth", "contig_percent_id",
                   "num_reads", "read_gene_id", "read_coverage_breadth", "read_coverage_depth"]]

        final_df.sort_index(inplace=True)
        final_df.dropna(subset=['drug_class'], inplace=True)
        final_df.to_csv("primary_AMR_report.tsv", sep='\t', index_label="gene_name")
        CODE
    >>>

    output {
        File contigs_tsv = "main_output.tsv"
        File reads_tsv = "kma_output.tsv"
        File final_summary_tsv = "comprehensive_AMR_metrics.tsv"
        File bigtable_tsv = "bigtable_report.tsv"
        File synthesized_report_tsv = "primary_AMR_report.tsv"
    }
     runtime {
        docker: "~{docker}"
    }
}


task motus {

    input {
        File fastq_1
        File fastq_2
        Int threads
        String docker
        String sample_name
    }

    command <<<
        set -euxo pipefail

        motus profile \
        -f ~{fastq_1} \
        -r ~{fastq_2} \
        -n ~{sample_name} \
        -t ~{threads} \
        -o motus.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File out = "motus.txt"
    }
}

task abricate {

    input {
        File contigs
        String docker
    }

    command <<<
        set -euxo pipefail

        abricate ~{contigs} > abricate.tsv
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File out = "abricate.tsv"
    }
}

task abritamr {

    input {
        File contigs
        String docker
    }

    command <<<
        set -euxo pipefail

        abritamr run --contigs ~{contigs}
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File? out = "./abritamr/abritamr.txt"
        File? matches = "./abritamr/summary_matches.txt"
        File? virulence = "./abritamr/summary_virulence.txt"
    }
}

task blast {

    input {
        File contigs
        Int threads
        String docker
    }

    command <<<
        set -ex -o pipefail
        touch blast_res.txt
        exit
        /opt/blast/ncbi-blast-2.9.0+/bin/blastn \
        -query ~{contigs} \
        -db /home/cromwell/mycob-ref/blast/nt_prok \
        -out blast_res.txt \
        -num_threads ~{threads} \
        -outfmt "6 qaccver saccver pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qlen slen qcovs stitle" \
        -max_target_seqs 10 \
        -evalue 1e-50
    >>>

    runtime {
        docker: "~{docker}"
        backend: "DockerAnnovar"
    }

    output {
        File blast_res = "blast_res.txt"
    }
}