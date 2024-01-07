version 1.0

task pangolin_one_sample {
    meta {
        description: "Pangolin classification of one SARS-CoV-2 sample."
    }
    input {
        File    genome_fasta
        Int?    min_length
        Float?  max_ambig
        Boolean inference_usher
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.3-pangolearn-2021-06-15"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.11-pangolearn-2021-08-09"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.14-pangolearn-2021-09-28"
        #String  docker = "staphb/pangolin:latest"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.14-pangolearn-2021-10-13"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.16-pangolearn-2021-11-18"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.16-pangolearn-2021-11-25"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.17-pangolearn-2022-01-05"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:3.1.20-pangolearn-2022-02-02"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.0.4-pangolearn-2022-03-22"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.0.5-pangolearn-2022-04-09"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.0.6-pangolearn-2022-05-23"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.1.2-pangolearn-2022-07-09"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.1.3-data-1.15.1"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.1.3-data-1.16"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.1.3-data-1.17"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.2-data-1.18.1.1"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.2-data-1.19"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.3-data-1.20"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.3-data-1.21"
        #String  docker = "cr.yandex/crp0pd4faec3e1qvam30/pangolin:4.3-data-1.22"
        String  docker
    }
    String basename = basename(genome_fasta, ".fasta")
    command <<<
        date | tee DATE
        conda list -n pangolin | grep "usher" | awk -F ' +' '{print$1, $2}' | tee VERSION_PANGO_USHER
        set -ex
        pangolin -v | tee VERSION_PANGOLIN
        pangolin -pv | tee VERSION_PANGOLEARN

        pangolin "~{genome_fasta}" \
            ~{true='--analysis-mode usher' false='--analysis-mode pangolearn' inference_usher} \
            --outfile "~{basename}.pangolin_report.csv" \
            ~{"--min-length " + min_length} \
            ~{"--max-ambig " + max_ambig} \
            --alignment \
            --verbose

        #cp sequences.aln.fasta "~{basename}.pangolin_msa.fasta"
        cp alignment.fasta "~{basename}.pangolin_msa.fasta"
        python3 <<CODE
        import csv
        #grab output values by column header
        with open("~{basename}.pangolin_report.csv", 'rt') as csv_file:
            for line in csv.DictReader(csv_file):
                with open("VERSION", 'wt') as outf:
                    pangolin_version=line["pangolin_version"]
                    version=line["version"]
                    outf.write(f"pangolin {pangolin_version}; {version}")
                with open("PANGO_LINEAGE", 'wt') as outf:
                    outf.write(line["lineage"])
                with open("PANGOLIN_CONFLICTS", 'wt') as outf:
                    outf.write(line["conflict"])
                with open("SCORPIO_CALL", 'wt') as outf:
                    outf.write(line["scorpio_call"])
                with open("PANGOLIN_NOTES", 'wt') as outf:
                    outf.write(line["note"])
                break
        CODE
    >>>
    runtime {
        docker: "~{docker}"
    }
    output {
        String     date                   = read_string("DATE")
        String     version                = read_string("VERSION")
        String     pango_lineage          = read_string("PANGO_LINEAGE")
        String     scorpio_call           = read_string("SCORPIO_CALL")
        String     pangolin_conflicts     = read_string("PANGOLIN_CONFLICTS")
        String     pangolin_notes         = read_string("PANGOLIN_NOTES")
        String     pangolin_usher_version = read_string("VERSION_PANGO_USHER")
        String     pangolin_version       = read_string("VERSION_PANGOLIN")
        String     pangolearn_version     = read_string("VERSION_PANGOLEARN")
        File       pango_lineage_report   = "${basename}.pangolin_report.csv"
        File       msa_fasta              = "~{basename}.pangolin_msa.fasta"
    }
}