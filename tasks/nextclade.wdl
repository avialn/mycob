version 1.0

task Nextclade {

    input {
        Array[File] ha_fasta
        Array[File] na_fasta
        String ref_HA="HA_wis_67_2022"
        String ref_NA="NA_wis_67_2022"
        String docker
    }

    Int len_ha = length(ha_fasta)
    Int len_na = length(na_fasta)

    command <<<
        set -ex -o pipefail

        if [ ~{len_ha} == 0 ]; then
            echo "No fasta file for HA segment in irma results!"
        else
            /opt/nextclade/nextclade run \
            --input-tree=/opt/nextclade/tree.json \
            --input-gene-map=/opt/nextclade/~{ref_HA}.gff \
            --input-qc-config=/opt/nextclade/qc.json \
            --input-pcr-primers=/opt/nextclade/primers.csv \
            --input-virus-properties=/opt/nextclade/virus_properties.json \
            --input-ref=/opt/nextclade/~{ref_HA}.fasta \
            --output-tsv HA_nextclade.tsv \
            ~{sep=' ' ha_fasta}
        fi

        if [ ~{len_na} == 0 ]; then
            echo "No fasta file for NA segment in irma results!"
        else
            /opt/nextclade/nextclade run \
            --input-tree=/opt/nextclade/tree.json \
            --input-gene-map=/opt/nextclade/~{ref_NA}.gff \
            --input-qc-config=/opt/nextclade/qc.json \
            --input-pcr-primers=/opt/nextclade/primers.csv \
            --input-virus-properties=/opt/nextclade/virus_properties.json \
            --input-ref=/opt/nextclade/~{ref_NA}.fasta \
            --output-tsv NA_nextclade.tsv \
            ~{sep=' ' na_fasta}
        fi
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File? HA_nextclade = "HA_nextclade.tsv"
        File? NA_nextclade = "NA_nextclade.tsv"
    }
}


task NextcladeParse {

    input {
        File? HA_tsv
        File? NA_tsv
        File report
        File antigenic_frame
        String docker
    }

    command <<<
        set -e
        python3 <<CODE

        import pandas as pd
        import json
        from typing import Dict, Optional, Tuple
        import logging
        import os.path
        from pathlib import Path
        import re

        nextclade_HA="~{HA_tsv}"
        nextclade_NA="~{NA_tsv}"
        report="~{report}"
        user_frame="~{antigenic_frame}"

        default_frame=[128, 129, 156, 157, 158, 159, 160, 162, 163, 164, 165, \
               166, 167, 187, 188, 189, 190, 191, 192, 193, 194, 195, \
               196, 197, 198, 169, 170, 171, 172, 173, 206, 207, 208, \
               238, 239, 240, 141, 142, 143, 144, 145, 224, 225, 74, 75, 76, 77, 78, 79]

        if not os.path.isfile(user_frame):
            antigenic_frame=default_frame
        else:
            read_frame = [i for i in Path(user_frame).read_text(encoding="utf-8").replace("\n", " ").split()]
            antigenic_frame = [int(item) for item in read_frame]

        def write_report_HA (
            summary: Optional[Dict],
            report_file: Optional[str]
        ):
            with open(report_file) as json_file:
                data = json.load(json_file)
            data["aa_substitutions_in_HA"]=summary["total_mut"]
            data["aa_substitutions_in_antigenic_sites_of_HA"]=summary["ant_mut"]
            with open(report_file, 'w') as outfile:
                json.dump(data, outfile, indent=4)

        def report_HA (
            nextclade_results: Optional[str],
            report_file: Optional[str],
            antigenic_frame=default_frame
        ):
            if not os.path.isfile (nextclade_results):
                with open(report_file) as json_file:
                    data = json.load(json_file)
                data["aa_substitutions_in_HA"] = "-"
                data["aa_substitutions_in_antigenic_sites_of_HA"] = "-"
                with open(report_file, 'w') as outfile:
                    json.dump(data, outfile, indent=4)
            else:
                antigenic_frame.sort()
                df = pd.read_csv(nextclade_results, sep='\t')
                if df["aaSubstitutions"].isna().item():
                    parse = "-"
                    ant_mut = "-"
                else:
                    parse = df["aaSubstitutions"].str.split(',').item()
                    ant_mut = []
                    for mut in parse:
                        dig = int(re.search( r'(?<=:(\w|\*))[0-9]*', mut).group())
                        if dig in antigenic_frame:
                            ant_mut.append(mut)
                    parse = " ".join(parse)
                    if len(ant_mut) == 0:
                        ant_mut = "-"
                    else:
                        ant_mut = " ".join(ant_mut)
                results = {}
                results["total_mut"] = parse
                results["ant_mut"] = ant_mut
                write_report_HA(results, report_file)

        def report_NA (
            nextclade_results: Optional[str] ,
            report_file: Optional[str]
        ):
            if not os.path.isfile(nextclade_results):
                with open(report_file) as json_file:
                    data = json.load(json_file)
                data["aa_substitutions_in_NA"] = "-"
                with open(report_file, 'w') as outfile:
                    json.dump(data, outfile, indent=4)
            else:
                df = pd.read_csv(nextclade_results, sep='\t')
                if df["aaSubstitutions"].isna().item():
                    logging.error(f"{df['errors'].item()}")
                    parse = "-"
                else:
                    parse = df["aaSubstitutions"].str.split(',').item()
                    parse = " ".join(parse)
                with open(report_file) as json_file:
                    data = json.load(json_file)
                data["aa_substitutions_in_NA"] = parse
                with open(report_file, 'w') as outfile:
                    json.dump(data, outfile, indent=4)

        report_HA (nextclade_HA, report, antigenic_frame)

        report_NA (nextclade_NA, report)

        CODE
    >>>

    runtime {
        docker:"~{docker}"
    }
}
