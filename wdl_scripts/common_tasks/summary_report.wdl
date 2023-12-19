version 1.0

task SummaryReport {

    input {
        String sample_name
        File? irma_report
        File? snpeff_HA_report
        File? snpeff_NA_report
        File? nextclade_HA_report
        File? nextclade_NA_report
        String docker
    }

    command <<<
        set -ex -o pipefail

        python3 <<CODE

        import json

        def merge_json_files(*json_files):
            merged_json = {}
            for json_file in json_files:
                if json_file:
                    with open(json_file, "r") as file:
                        data = json.load(file)
                        merged_json.update(data)
            return merged_json


        merged_json = merge_json_files ("~{irma_report}", "~{snpeff_HA_report}", "~{snpeff_NA_report}", "~{nextclade_HA_report}", "~{nextclade_NA_report}")

        with open("~{sample_name}_summary_report.json", "w") as outfile:
            json.dump(merged_json, outfile, indent=4)

        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report_json = "~{sample_name}_summary_report.json"
    }
}

