version 1.0

task Nextclade {

    input {
        Array[File] fasta
        String ref_name
        String docker
    }

    Int len_array = length(fasta)

    command <<<
        set -ex -o pipefail

        if [ ~{len_array} != 0 ]; then
          proceed="yes"
          /opt/nextclade/nextclade run \
          --input-tree=/opt/nextclade/tree.json \
          --input-gene-map=/opt/nextclade/~{ref_name}.gff \
          --input-qc-config=/opt/nextclade/qc.json \
          --input-pcr-primers=/opt/nextclade/primers.csv \
          --input-virus-properties=/opt/nextclade/virus_properties.json \
          --input-ref=/opt/nextclade/~{ref_name}.fasta \
          --output-tsv nextclade.tsv \
          ~{sep=' ' fasta}
          cov=$(awk -F'\t' 'NR == 2 {print $19}' nextclade.tsv)
        fi

        echo $proceed > proceed.txt
        echo $cov > coverage_percentage.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File? nextclade_tsv = "nextclade.tsv"
        String proceed = read_string("proceed.txt")
        String coverage_percentage = read_string("coverage_percentage.txt")
    }
}


task NextcladeParse {

    input {
        File? nextclade_tsv
        File? antigenic_frame
        String search_antigenic_mut
        String ref_name
        String docker
    }

    command <<<
        set -ex -o pipefail

        python3 <<CODE

        import pandas as pd
        import json
        from typing import Dict, Optional, Tuple, List
        import os.path
        from pathlib import Path
        import re

        default_frame=[128, 129, 156, 157, 158, 159, 160, 162, 163, 164, 165, \
               166, 167, 187, 188, 189, 190, 191, 192, 193, 194, 195, \
               196, 197, 198, 169, 170, 171, 172, 173, 206, 207, 208, \
               238, 239, 240, 141, 142, 143, 144, 145, 224, 225, 74, 75, 76, 77, 78, 79]


        if not os.path.isfile("~{antigenic_frame}"):
            antigenic_frame=default_frame
        else:
            read_frame = [i for i in Path("~{antigenic_frame}").read_text(encoding="utf-8").replace("\n", " ").split()]
            antigenic_frame = [int(item) for item in read_frame]


        def parse_mutations (
            nextclade_results: Optional[str],
            antigenic_frame: Optional[List[int]]
        ):
            antigenic_frame.sort()
            df = pd.read_csv(nextclade_results, sep='\t')
            if df["aaSubstitutions"].isna().item():
                mut = "-"
                ant_mut = "-"
            else:
                mut = df["aaSubstitutions"].str.split(',').item()
                ant_mut = []
                for m in mut:
                    dig = int(re.search( r'(?<=:(\w|\*))[0-9]*', m).group())
                    if dig in antigenic_frame:
                        ant_mut.append(m)
                mut = " ".join(mut)
                if len(ant_mut) == 0:
                    ant_mut = "-"
                else:
                    ant_mut = " ".join(ant_mut)
            return mut, ant_mut


        def write_report (
            ref_name: Optional[str],
            nextclade_results: Optional[str],
            search_antigenic_mut: Optional[str],
            antigenic_frame: Optional[List[str]]

        ):
            mut, ant_mut = parse_mutations (nextclade_results, antigenic_frame)
            if search_antigenic_mut == "yes":
                with open("nextclade.json", 'w') as json_file:
                    data = {
                        f"aaSubstitutions_to_{ref_name}_nextclade" : mut,
                        f"antigenic_site_aaSubstitutions_to_{ref_name}_nextclade" : ant_mut
                    }
                    json.dump(data, json_file, indent=4)
            else:
                with open("nextclade.json", 'w') as json_file:
                    data = {
                        f"aaSubstitutions_to_{ref_name}_nextclade" : mut
                    }
                    json.dump(data, json_file, indent=4)

        write_report ("~{ref_name}", "~{nextclade_tsv}", "~{search_antigenic_mut}", antigenic_frame)

        CODE
    >>>

    runtime {
        docker:"~{docker}"
    }

    output {
       File nextclade_json = "nextclade.json"
    }
}
