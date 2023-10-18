version 1.0

task Irma {

    input {
        File trim_R1
        File trim_R2
        String module_irma
        String docker
    }

    command <<<
        set -ex -o pipefail

        /opt/irma/flu-amd/IRMA ~{module_irma} ~{trim_R1} ~{trim_R2} irma_res

        # determine if assemly was successful
        if compgen -G "irma_res/*.fasta"; then
            proceed="yes"
        else
            proceed="no IRMA assembly generated"
        fi

        echo $proceed > proceed.txt

        if [ $proceed == "yes" ] && [ ~{module_irma} == "FLU" ]; then
            TYPE=$(basename $(find irma_res/*.fasta | head -n 1 ) | cut -d "_" -f 1)

            ## Type == A then grab subtype if exists
            if [ $TYPE == "A" ]; then
                if compgen -G "irma_res/*_HA*.fasta"; then
                    HA_SUBTYPE=$(basename $(find irma_res/*HA*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    HA_SUBTYPE="-"
                fi

                if compgen -G "irma_res/*_NA*.fasta"; then
                    NA_SUBTYPE=$(basename $(find irma_res/*NA*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    NA_SUBTYPE="-"
                fi

            else
                HA_SUBTYPE="-"
                NA_SUBTYPE="-"
            fi
        fi

        echo ${TYPE}${HA_SUBTYPE}${NA_SUBTYPE} > irma_flu_type.txt

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        String proceed = read_string("proceed.txt")
        String flu_type = read_string("irma_flu_type.txt")
        Array[File] irma_fasta  = glob("irma_res/*.fasta")
        Array[File] ha_fasta = select_first([glob("irma_res/*HA*.fasta"), " "])
        Array[File] na_fasta = select_first([glob("irma_res/*NA*.fasta"), " "])
        File? irma_qc = "irma_res/tables/READ_COUNTS.txt"
    }
}

task IrmaQC {

    input {
        String sample_name
        String module_irma
        File? irma_qc
        String docker
        File reference_fasta
        Array[File] irma_fasta
    }

    command <<<
        set -ex -o pipefail

        # removing \n in the reference fasta
        awk '/^>/ {if (NR!=1) {printf("\n%s\n", $0)} \
        else {printf("%s\n",$0)}; next} { printf("%s",$0)} END {printf("\n")}' \
        ~{reference_fasta} > ref.fasta

        echo "fasta_name	ref_length    fasta_length  length_assembling_ratio" > ~{sample_name}_~{module_irma}_irma_length.txt

        # finding the ratio of the assembled fasta length to the reference fasta length
        declare -A data
        for fasta_file in ~{sep=" " irma_fasta}; do
          fasta_name=$(basename "${fasta_file}" .fasta)
          ref_length=$(grep -A 1 "^>${fasta_name}$" ref.fasta \
          | awk '/^[^>]/ {seq = seq $0} END {print length(seq)}')
          fasta_length=$(cat $fasta_file | awk '/^[^>]/ {seq = seq $0} END {print length(seq)}')
          ratio=$(printf "%.2f%%\n" $(echo "$fasta_length / $ref_length * 100" | bc -l))
          data[$fasta_name]=$ratio
          echo $fasta_name $ref_length $fasta_length $ratio >> ~{sample_name}_~{module_irma}_irma_length.txt
        done

        # writing json short dictionary to zip it in the irma_qc report
        json_data="{"
        for key in "${!data[@]}"; do
          json_data+="\"$key\": \"${data[$key]}\", "
        done
        json_data="${json_data%, *} }"

        echo $json_data > irma_length_ratio.json

        python3 <<CODE

        import pandas as pd
        import json

        data = pd.read_csv ("~{irma_qc}", sep = "\t")

        reads_total = data.loc [data["Record"] == "1-initial"]["Reads"].item()
        reads_passQC = data.loc [data["Record"] == "2-passQC"]["Reads"].item()
        reads_match = data.loc [data["Record"] == "3-match"]["Reads"].item()

        data_fasta = data[~data["Record"].str.startswith(('0', '1', '2', '3'))]
        data_fasta["Record"] = data_fasta["Record"].str[2:]

        # finding the virus type to which the majority of reads align
        sorted_data = data_fasta.sort_values(by='Reads', ascending=False)
        type_irma = sorted_data.iloc[0]["Record"]

        with open ("~{sample_name}_~{module_irma}_irma_type.txt", 'w') as file:
            file.write (type_irma )

        qc_dict = dict(zip(data_fasta["Record"], data_fasta["Reads"]))

        with open("irma_length_ratio.json", "r") as ratio_file:
            irma_length_ratio = json.load(ratio_file)

        with open("~{sample_name}_~{module_irma}_irma_qc.json", 'w') as f:
            data = {
                "sample_name": "~{sample_name}",
                "total_reads": reads_total,
                "reads_passed_QC": reads_passQC,
                "reads_matched": reads_match,
                "fasta_assembled": qc_dict,
                "length_assembling_ratio": irma_length_ratio
            }
            json.dump(data, f, indent=4)
        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File irma_qc_json = "~{sample_name}_~{module_irma}_irma_qc.json"
        File irma_length_txt = "~{sample_name}_~{module_irma}_irma_length.txt"
        String irma_type = read_string("~{sample_name}_~{module_irma}_irma_type.txt")
    }
}

task Blast {

    input {
        String sample_name
        String module_irma
        Array[File] fasta_files
        String db
        Int num_alignments=100
        Int threads
        String docker
    }

    command <<<
        set -ex -o pipefail

        cat ~{sep=" " fasta_files} > tmp.fasta
        sed 's/\./N/g' tmp.fasta > big.fasta

        /opt/blast/ncbi-blast-2.9.0+/bin/blastn \
        -query big.fasta \
        -db /opt/blast/~{db} \
        -out ~{sample_name}_~{module_irma}_blast_res \
        -num_threads ~{threads} \
        -outfmt "6 qaccver saccver pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qlen slen qcovs stitle" \
        -num_alignments ~{num_alignments} \
        -evalue 1e-6

        rm tmp.fasta big.fasta
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File blast_res = "~{sample_name}_~{module_irma}_blast_res"
    }
}

task Report {

    input {
        String sample_name
        String module_irma
        File? blast_res
        String irma_type
        Float pi_thresh=0.8
        Int min_aln_len=50
        String docker
    }

    command <<<
        set -ex -o pipefail
        python3 <<CODE

        import pandas as pd
        import numpy as np
        from typing import Dict, List, Optional, Tuple
        import json
        import re
        pd.options.mode.chained_assignment = None

        blast_cols = [
            ("qaccver", "category"),
            ("saccver", str),
            ("pident", float),
            ("length", "uint16"),
            ("mismatch", "uint16"),
            ("gapopen", "uint16"),
            ("qstart", "uint16"),
            ("qend", "uint16"),
            ("sstart", "uint16"),
            ("send", "uint16"),
            ("evalue", np.float16),
            ("bitscore", np.float16),
            ("qlen", "uint16"),
            ("slen", "uint16"),
            ("qcovs", np.float16),
            ("stitle", str),
        ]

        md_cols = [
            ("accession", str),
            ("host", "category"),
            ("segment", "category"),
            ("subtype", "str"),
            ("country", "category"),
            ("date", "category"),
            ("seq_length", "uint16"),
            ("virus_name", "category"),
            ("age", "category"),
            ("gender", "category"),
            ("group_id", "category"),
        ]

        def write_irma_report(
            report_file: Optional[str],
            sample_name: Optional[str],
            module_irma: Optional[str],
            irma_type: Optional[str]
        ):
            irma_report = {
                "sample_name": sample_name,
                "module_irma": module_irma,
                "irma_virus_type": irma_type
            }
            with open(report_file, "w") as json_file:
                json.dump(irma_report, json_file, indent=4)

        write_irma_report ("~{sample_name}_~{module_irma}_report.json", "~{sample_name}", "~{module_irma}", "~{irma_type}")

        # blast virus name (excluding cases of COVID and influenza)
        def parse_blast (
            blast_results: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int]
        ) -> Optional[Dict]:
            df = pd.read_csv(
                blast_results, sep='\t',
                names=[name for name, coltype in blast_cols],
                dtype={name: coltype for name, coltype in blast_cols},
            )
            if df.empty:
                results_summary = {
                    "blast_virus_name": "Blast_results_are_empty"
                }
            else:
                df['stitle'] = df['stitle'].str.replace('|', '', regex=True)
                df_filtered = df.loc[
                    (df["pident"] >= (pident_threshold * 100)) & (df["length"] >= min_aln_length)
                ]
                if df_filtered.empty:
                    results_summary = {
                        "blast_virus_name": "Try to modify pident_threshold or min_aln_length"
                    }
                else:
                    df_filtered.loc[:,"qaccver"] = pd.Categorical(df_filtered.loc[:,"qaccver"])
                    df_filtered = df_filtered.sort_values(by="bitscore", ascending= False).reset_index()
                    full_virus_name=df_filtered.loc [0, "stitle"]
                    results_summary = {
                        "blast_virus_name": full_virus_name
                    }
            return results_summary

        def write_blast_report (
            report_file: Optional[str],
            blast_results: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int]
        ):
            blast_report = parse_blast (blast_results, pident_threshold, min_aln_length)
            with open(report_file) as file:
                data = json.load(file)
            data["blast_virus_name"] = blast_report["blast_virus_name"]
            with open(report_file, 'w') as outfile:
                json.dump(data, outfile, indent=4)

        ################ FLU
        def parse_blast_flu (
            blast_results: Optional[str],
            flu_metadata: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int]
        ) -> Optional[Dict]:
            df_md = pd.read_csv(
            flu_metadata,
            sep="\t",
            names=[name for name, _ in md_cols],
            dtype={name: t for name, t in md_cols},
            )
            unique_subtypes = df_md.subtype.unique()
            unique_subtypes = unique_subtypes[~pd.isna(unique_subtypes)]
            regex_subtype_pattern = r"\((H\d+N\d+|" + "|".join(list(unique_subtypes)) + r")\)"
            df = pd.read_csv(
                blast_results, sep='\t',
                names=[name for name, coltype in blast_cols],
                dtype={name: coltype for name, coltype in blast_cols},
            )
            if df.empty:
                results_summary = {
                    "blast_virus_name": "Blast_results_are_empty"
                }
            else:
                df_filtered = df.loc[
                    (df["pident"] >= (pident_threshold * 100)) & (df["length"] >= min_aln_length)
                ]
                if df_filtered.empty:
                    results_summary = {
                        "blast_virus_name": "Try to modify pident_threshold or min_aln_length"
                    }
                else:
                    df_filtered["accession"] = df_filtered["saccver"].str.strip()
                    df_filtered["qaccver"] = pd.Categorical(df_filtered["qaccver"])
                    df_filtered.loc[:, "subtype_from_match_title"] = (
                    df_filtered.loc[:, "stitle"].str.extract(regex_subtype_pattern).astype("category")
                    )
                    segments = df_filtered.qaccver.unique()
                    influenza_segment = {}
                    for segment in segments:
                        if "_PB2" in segment:
                            influenza_segment[segment] = 1
                        if "_PB1" in segment:
                            influenza_segment[segment] = 2
                        if "_PA" in segment:
                            influenza_segment[segment] = 3
                        if "_HA" in segment:
                            influenza_segment[segment] = 4
                        if "_NP" in segment:
                            influenza_segment[segment] = 5
                        if "_NA" in segment:
                            influenza_segment[segment] = 6
                        if "_M" in segment:
                            influenza_segment[segment] = 7
                        if "_NS" in segment:
                            influenza_segment[segment] = 8
                    df_filtered["qaccver"] = df_filtered["qaccver"].map(influenza_segment)
                    df_filtered = df_filtered.sort_values(
                        by=["qaccver", "bitscore"], ascending=[True, False]
                    ).set_index("qaccver")
                    reg_pat_H = r"(H\d)"
                    reg_pat_N = r"(N\d)"
                    subtype_HA, virus_name_HA = parse_HA_NA(4, df_filtered, df_md, reg_pat_H)
                    subtype_NA, virus_name_NA = parse_HA_NA(6, df_filtered, df_md, reg_pat_N)
                    subtype = get_subtype_value(subtype_HA, subtype_NA)
                    results_summary = {
                        "blast_subtype": subtype,
                        "virus_name_HA": virus_name_HA,
                        "virus_name_NA": virus_name_NA
                    }
            return results_summary

        def parse_HA_NA (
            seg: Optional[int],
            df: Optional[pd.DataFrame],
            df_md: Optional[pd.DataFrame],
            reg_pat: Optional[str]
        ) -> Optional[Tuple]:
            if seg in df.index:
                df_seg = df.loc[seg, :]
                tmp = df_seg["subtype_from_match_title"].value_counts()
                if len(tmp) == 0:
                    subtype = None
                else:
                    result = re.search(reg_pat, tmp.index[0])
                    subtype = result.group(0)
                df_merge = pd.merge(df_seg, df_md, on="accession", how="left")
                if df_merge.empty:
                    virus_name = "-"
                else:
                    virus_name = df_merge.loc[0, "virus_name"]
                return subtype, virus_name
            else:
                subtype = "-"
                virus_name = "-"
                return subtype, virus_name

        def get_subtype_value (
            subtype_HA: Optional[str],
            subtype_NA: Optional[str]
        ) -> Optional[str]:
            subtype = ""
            if subtype_HA is None and subtype_NA is None:
                subtype = "-"
            elif subtype_HA is not None and subtype_NA is None:
                H: str = subtype_HA
                subtype = f"{H}" if H != "" else "-"
            elif subtype_HA is None and subtype_NA is not None:
                N: str = subtype_NA
                subtype = f"{N}" if N != "" else "-"
            else:
                H: str = subtype_HA
                N: str = subtype_NA
                if H == "" and N == "":
                    subtype = "-"
                else:
                    if H != "":
                        H = f"{H}"
                    if N != "":
                        N = f"{N}"
                    subtype = f"{H}{N}"
            return subtype

        def write_flu_blast_report (
            report_file: Optional[str],
            blast_results: Optional[str],
            flu_metadata: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int]
        ):
            blast_report = parse_blast_flu (
            blast_results,
            flu_metadata,
            pident_threshold,
            min_aln_length
            )
            with open(report_file) as file:
                data = json.load(file)
            data["blast_type"] = blast_report ["blast_subtype"]
            data["blast_virus_name_HA"] = blast_report ["virus_name_HA"]
            data["blast_virus_name_NA"] = blast_report ["virus_name_NA"]
            with open(report_file, 'w') as outfile:
                json.dump(data, outfile, indent=4)

        if ("~{module_irma}" != "CoV") and ("~{module_irma}" != "FLU"):
            write_blast_report ("~{sample_name}_~{module_irma}_report.json", "~{blast_res}", ~{pi_thresh}, ~{min_aln_len})
        elif ("~{module_irma}" == "FLU"):
            flu_metadata = "/genomeset.dat"
            write_flu_blast_report ("~{sample_name}_~{module_irma}_report.json", "~{blast_res}", flu_metadata,  ~{pi_thresh}, ~{min_aln_len})

        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report = "~{sample_name}_~{module_irma}_report.json"
    }
}
