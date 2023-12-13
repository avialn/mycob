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

task IrmaGut {

    input {
        File trim_R1
        File trim_R2
        String module_irma
        String docker
    }

    command <<<
        set -ex -o pipefail

        /opt/flu-amd/IRMA ~{module_irma} ~{trim_R1} ~{trim_R2} irma_res

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

        if [ $proceed == "yes" ] && [ ~{module_irma} == "ROTAVIRUS" ]; then
            TYPE=$(basename $(find irma_res/*.fasta | head -n 1 ) | cut -d "_" -f 1)

            ## Type == A then grab subtype if exists
            if [ $TYPE == "A" ]; then
                if compgen -G "irma_res/*_VP7*.fasta"; then
                    G_SUBTYPE=$(basename $(find irma_res/*VP7*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    G_SUBTYPE="-"
                fi

                if compgen -G "irma_res/*_VP4*.fasta"; then
                    P_SUBTYPE=$(basename $(find irma_res/*VP4*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    P_SUBTYPE="-"
                fi

            else
                G_SUBTYPE="-"
                P_SUBTYPE="-"
            fi
        fi

        echo ${TYPE}${HA_SUBTYPE}${NA_SUBTYPE} > irma_flu_type.txt
        echo ${TYPE}${G_SUBTYPE}${P_SUBTYPE} > irma_rotavirus_type.txt

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        String proceed = read_string("proceed.txt")
        String flu_type = read_string("irma_flu_type.txt")
        String rotavirus_type = read_string("irma_rotavirus_type.txt")
        Array[File] irma_fasta  = glob("irma_res/*.fasta")
        Array[File] ha_fasta = select_first([glob("irma_res/*HA*.fasta"), " "])
        Array[File] na_fasta = select_first([glob("irma_res/*NA*.fasta"), " "])
        Array[File] g_fasta = select_first([glob("irma_res/*_VP7*.fasta"), " "])
        Array[File] p_fasta = select_first([glob("irma_res/*_VP4*.fasta"), " "])
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

        VERSION=$(ls /opt/blast/ | grep "ncbi")

        /opt/blast/$VERSION/bin/blastn \
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
        File? rotavirus_metadata
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

        # parse virus type by irma fasta (for flu, rotavirus and all others)
        def write_irma_report(
            report_file: Optional[str],
            sample_name: Optional[str],
            module_irma: Optional[str],
            irma_type: Optional[str]
        ):
            if module_irma == "FLU":
                A_or_B = irma_type[0]
                if len(irma_type) > 1:
                    subtype = irma_type[1:]
                else:
                    subtype = "-"
                if A_or_B == "A":
                    irma_report = {
                    "sample_name": sample_name,
                    "module_irma": module_irma,
                    "irma_virus_type": f"Influenza {A_or_B} {subtype}"
                    }
                else:
                    irma_report = {
                    "sample_name": sample_name,
                    "module_irma": module_irma,
                    "irma_virus_type": f"Influenza {A_or_B}"
                    }
            elif module_irma == "ROTAVIRUS":
                type_virus = irma_type[0]
                if len(irma_type) > 1:
                    subtype = irma_type[1:]
                    subtype = re.sub(r'(P\d+)', r'[\1]', subtype)
                else:
                    subtype = "-"
                if type_virus == "A":
                    irma_report = {
                    "sample_name": sample_name,
                    "module_irma": module_irma,
                    "irma_virus_type": f"Rotavirus {type_virus} {subtype}"
                    }
                else:
                    irma_report = {
                    "sample_name": sample_name,
                    "module_irma": module_irma,
                    "irma_virus_type": f"Rotavirus {type_virus}"
                    }
            else:
                virus_name = irma_type.replace("_", " ")
                irma_report = {
                    "sample_name": sample_name,
                    "module_irma": module_irma,
                    "irma_virus_type": virus_name
            }
            with open(report_file, "w") as json_file:
                json.dump(irma_report, json_file, indent=4)

        # parse virus name by blast results (NOT for covid or influenza or rotavirus)
        def parse_blast (
            blast_results: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int],
            irma_type: Optional[str]
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
                df_filtered = df[df.qaccver == irma_type].reset_index(drop=True)
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
                    full_virus_name=df_filtered.loc [0, "stitle"].replace('|', '').replace(', complete genome', '')
                    results_summary = {
                        "blast_virus_name": full_virus_name
                    }
            return results_summary

        def write_blast_report (
            report_file: Optional[str],
            blast_results: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int],
            irma_type: Optional[str]
        ):
            blast_report = parse_blast (blast_results, pident_threshold, min_aln_length, irma_type)
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
            min_aln_length: Optional[int],
            irma_type: Optional[str]
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
                    A_or_B = irma_type[0]
                    df_filtered["accession"] = df_filtered["saccver"].str.strip()
                    df_filtered["qaccver"] = pd.Categorical(df_filtered["qaccver"])
                    df_filtered.loc[:, "subtype_from_match_title"] = (
                    df_filtered.loc[:, "stitle"].str.extract(regex_subtype_pattern).astype("category")
                        )
                    segments = df_filtered.qaccver.unique()
                    df_merge = pd.merge(df_filtered, df_md, on="accession", how="left")

                    # find the names with the highest bitscore score for each segment
                    if df_merge.empty:
                        virus_name_allsegments = "-"
                    else:
                        names = []
                        for seg in segments:
                            df_sort_dropNaN = df_merge[df_merge.qaccver == seg]\
                            .sort_values(by=["bitscore"], ascending=[False])\
                            .dropna(subset=["virus_name"])\
                            .reset_index(drop=True)
                            if not df_sort_dropNaN.empty:
                                best_name = df_sort_dropNaN["virus_name"].iloc[0]
                                names.append(best_name)
                        # find the most common name between all segments
                        if len(names)>0:
                            virus_name_allsegments = max(names, key=names.count)
                        else:
                            virus_name_allsegments = "-"

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

                    if A_or_B == "A":
                        #subtype by all segments
                        subtype_counts = df_filtered['subtype_from_match_title'].value_counts()
                        if len(subtype_counts) == 0:
                            subtype_allsegments = "-"
                        else:
                            subtype_allsegments = subtype_counts.idxmax()

                        results_summary = {
                            "blast_subtype_HANA": subtype,
                            "blast_subtype_allsegments": subtype_allsegments,
                            "virus_name_HA": virus_name_HA,
                            "virus_name_NA": virus_name_NA,
                            "virus_name_allsegments": virus_name_allsegments
                        }

                    if A_or_B == "B":
                        results_summary = {
                            "virus_name_HA": virus_name_HA,
                            "virus_name_NA": virus_name_NA,
                            "virus_name_allsegments": virus_name_allsegments
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
                    df_dropNaN = df_merge.dropna(subset=["virus_name"]).reset_index(drop=True)
                    if not df_dropNaN.empty:
                        virus_name = df_dropNaN.loc[0, "virus_name"]
                    else:
                        virus_name = "-"
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
            min_aln_length: Optional[int],
            irma_type: Optional[str]
        ):
            with open(report_file) as file:
                data = json.load(file)
            blast_data = parse_blast_flu (blast_results, flu_metadata, pident_threshold, min_aln_length, irma_type)
            data.update(blast_data)
            with open(report_file, 'w') as outfile:
                json.dump(data, outfile, indent=4)

        ############### ROTAVIRUS

        def find_G_genotype(string):
            matches = re.findall(reg_pat_G, str(string))
            return ",".join(matches) if matches else None

        def find_P_genotype(string):
            matches = re.findall(reg_pat_P, str(string))
            return ",".join(matches) if matches else None

        def parse_blast_rotavirus (
            blast_results: Optional[str],
            rotavirus_metadata: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int],
            irma_type: Optional[str]
        ) -> Optional[Dict]:
            df_md = pd.read_csv(rotavirus_metadata)
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
                    type_virus = irma_type[0]
                    df_filtered["GenBank Accessions"] = df_filtered['saccver'].str.split('.').str[0]
                    df_filtered["qaccver"] = pd.Categorical(df_filtered["qaccver"])
                    segments = df_filtered.qaccver.unique()
                    df_merge = pd.merge(df_filtered, df_md, on="GenBank Accessions", how="left")

                   # find the names with the highest bitscore score for each segment
                    if df_merge.empty:
                        virus_name_allsegments = "-"
                    else:
                        names = []
                        for seg in segments:
                            df_seg_sort = df_merge[df_merge.qaccver == seg]\
                            .sort_values(by=["bitscore"], ascending=[False])\
                            .dropna(subset=["Genome Name"])\
                            .reset_index(drop=True)
                            if not df_seg_sort.empty:
                                best_name = df_seg_sort["Genome Name"].iloc[0]
                                names.append(best_name)
                        # find the most common name between all segments
                        if len(names)>0:
                            virus_name_allsegments = max(names, key=names.count)
                        else:
                            virus_name_allsegments = "-"

                    rotavirus_segment = {}
                    for segment in segments:
                        if "_VP1" in segment:
                            rotavirus_segment[segment] = 1
                        if "_VP2" in segment:
                            rotavirus_segment[segment] = 2
                        if "_VP3" in segment:
                            rotavirus_segment[segment] = 3
                        if "_VP4" in segment:
                            rotavirus_segment[segment] = 4
                        if "_NSP1" in segment:
                            rotavirus_segment[segment] = 5
                        if "_VP6" in segment:
                            rotavirus_segment[segment] = 6
                        if "_NSP3" in segment:
                            rotavirus_segment[segment] = 7
                        if "_NSP2" in segment:
                            rotavirus_segment[segment] = 8
                        if "_VP7" in segment:
                            rotavirus_segment[segment] = 9
                        if "_NSP4" in segment:
                            rotavirus_segment[segment] = 10
                        if "_NSP5" in segment:
                            rotavirus_segment[segment] = 11

                    df_merge["qaccver"] = df_merge["qaccver"].map(rotavirus_segment)
                    df_merge = df_merge.sort_values(
                        by=["qaccver", "bitscore"], ascending=[True, False]
                    ).set_index("qaccver")
                    reg_pat_G = r"(G\d{1,2})" #9 segment, _VP7
                    reg_pat_P = r"(P\[\d{1,2}\]|P\d{1,2}|P\(\d\))" #4 segment, _VP4
                    df_merge['G_genotype'] = df_merge['Genotype'].apply(find_G_genotype)
                    df_merge['P_genotype'] = df_merge['Genotype'].apply(find_P_genotype)
                    df_merge['P_genotype'] = df_merge['P_genotype'].str.replace(r'[\[\]]', '', regex=True)
                    subtype_G, virus_name_G = parse_G_P(9, df_merge, reg_pat_G)
                    subtype_P, virus_name_P = parse_G_P(4, df_merge, reg_pat_P)
                    subtype = get_subtype_value_GP(subtype_G, subtype_P)

                    if type_virus == "A":
                        # subtype by all segments
                        subtype_counts = df_merge['Genotype'].value_counts()
                        if len(subtype_counts) == 0:
                            subtype_allsegments = "-"
                        else:
                            subtype_allsegments = subtype_counts.idxmax()
                            subtype_allsegments = re.sub(r'(P\d+)', r'[\1]', subtype_allsegments)

                        results_summary = {
                            "blast_subtype_by_G[P]seg": subtype,
                            "blast_subtype_allsegments": subtype_allsegments,
                            "virus_name_G": virus_name_G,
                            "virus_name_[P]": virus_name_P,
                            "virus_name_allsegments": virus_name_allsegments
                        }

                    else:
                        results_summary = {
                            "virus_name_G": virus_name_G,
                            "virus_name_[P]": virus_name_P,
                            "virus_name_allsegments": virus_name_allsegments
                        }
            return results_summary

        def parse_G_P (
            seg: Optional[int],
            df: Optional[pd.DataFrame],
            reg_pat: Optional[str]
        ) -> Optional[Tuple]:
            if seg in df.index:
                df_seg = df.loc[seg, :]
                if seg == 9:
                    tmp = df_seg["G_genotype"].value_counts()
                else:
                    tmp = df_seg["P_genotype"].value_counts()
                if len(tmp) == 0:
                    subtype = None
                else:
                    subtype=tmp.index[0]
                if df.empty:
                    virus_name = "-"
                else:
                    df_dropNaN = df_seg.dropna(subset=["Genome Name"]).reset_index(drop=True)
                    if not df_dropNaN.empty:
                        virus_name = df_dropNaN.loc[0, "Genome Name"]
                    else:
                        virus_name = "-"
                return subtype, virus_name
            else:
                subtype = "-"
                virus_name = "-"
                return subtype, virus_name

        def get_subtype_value_GP (
            subtype_G: Optional[str],
            subtype_P: Optional[str]
        ) -> Optional[str]:
            subtype = ""
            if subtype_G is None and subtype_P is None:
                subtype = "-"
            elif subtype_G is not None and subtype_P is None:
                G: str = subtype_G
                subtype = f"{G}" if G != "" else "-"
            elif subtype_G is None and subtype_P is not None:
                P: str = subtype_P
                subtype = f"[{P}]" if P != "" else "-"
            else:
                G: str = subtype_G
                P: str = subtype_P
                if G == "" and P == "":
                    subtype = "-"
                else:
                    if G != "":
                        G = f"{G}"
                    if P != "":
                        P = f"{P}"
                    subtype = f"{G}{P}"
            return subtype

        def write_rotavirus_blast_report (
            report_file: Optional[str],
            blast_results: Optional[str],
            rotavirus_metadata: Optional[str],
            pident_threshold: Optional[float],
            min_aln_length: Optional[int],
            irma_type: Optional[str]
        ):
            with open(report_file) as file:
                data = json.load(file)
            blast_data = parse_blast_rotavirus (blast_results,
                             rotavirus_metadata,
                             pident_threshold,
                             min_aln_length,
                             irma_type
                            )
            data.update(blast_data)
            with open(report_file, 'w') as outfile:
                json.dump(data, outfile, indent=4)

        ####### WRITING
        # irma type
        write_irma_report ("~{sample_name}_~{module_irma}_report.json", \
                           "~{sample_name}", \
                           "~{module_irma}", \
                           "~{irma_type}" \
                           )

        # blast type
        if ("~{module_irma}" != "CoV") and ("~{module_irma}" != "FLU") and ("~{module_irma}" != "ROTAVIRUS"):
            write_blast_report ("~{sample_name}_~{module_irma}_report.json", \
                                "~{blast_res}", \
                                ~{pi_thresh}, \
                                ~{min_aln_len}, \
                                "~{irma_type}" \
                                )
        elif ("~{module_irma}" == "FLU"):
            flu_metadata = "/genomeset.dat"
            write_flu_blast_report ("~{sample_name}_~{module_irma}_report.json", \
                                   "~{blast_res}", \
                                   flu_metadata, \
                                   ~{pi_thresh}, \
                                   ~{min_aln_len}, \
                                   "~{irma_type}" \
                                   )
        elif ("~{module_irma}" == "ROTAVIRUS"):
            write_rotavirus_blast_report ("~{sample_name}_~{module_irma}_report.json", \
                                   "~{blast_res}", \
                                   "~{rotavirus_metadata}", \
                                   ~{pi_thresh}, \
                                   ~{min_aln_len}, \
                                   "~{irma_type}" \
                                   )

        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report = "~{sample_name}_~{module_irma}_report.json"
    }
}
