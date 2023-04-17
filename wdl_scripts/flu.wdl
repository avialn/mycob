version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run workflow.wdl -i inputs.json -o options.json

import "tasks.wdl" as tasks

workflow FluGenotyping {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Int threads_num
    }

    call tasks.TrimTask as trimmomatic_task {
        input:
            fastq_1=fastq_1,
            fastq_2=fastq_2,
            threads_num=threads_num
    }

    call tasks.IrmaTask as irma_task {
        input:
            trim_R1 = trimmomatic_task.trim_R1,
            trim_R2 = trimmomatic_task.trim_R2
    }

    # даже есди нет ни одного фаста файла собранного ирмой,
    # пустой массив  if (defined(IrmaTask.irma_fasta)) == true
    # а значит остановим воркфлоу перед запуском бласта по условию ~{len_array} == 0

    call tasks.BlastTask as blast_task {
        input:
            fasta_files=irma_task.irma_fasta,
            threads_num=threads_num
    }

    # если ирма выдала хоть один фаста файл -> массив фаст не был пуст а значит
    # BlastTask отработал и выдал blast_res
    # поэтому в любом случе надо вызывать  BlastParse

    call tasks.BlastParse as blast_parse_task {
        input:
            blast_row=blast_task.blast_res,
            sample_name=sample_name,
    }

    # передаем массивы, в каждом лежит по одной ссылке на фаста файл - ha или na сегмент
    # массивы получились на выходе функции select_first
    # функция select_first нужна была для передачи фаста файла который мог быть, а мог и не быть
    # при отсутвии фаста файла для NA сегмента команда
    # File? na_fasta = glob("irma_res/*NA*.fasta") выдвет следующую ошибку
    #Failed to evaluate job outputs:
    #Bad output 'IrmaTask.na_fasta': IllegalArgumentException:
    #No coercion defined from wom value(s) '[]' of type 'Array[Nothing]' to 'File?'.
    # поэтому обернули в select_first

    #Если в аутпут ирмы передать файл: File? ha_fasta = write_lines(glob("irma_res/*HA*.fasta"))
    # В инпуте некстклада обернуть его в ha_fasta = select_first([IrmaTask.ha_fasta])
    # нексклад считывает ha_fasta с ошибками и падает

    call tasks.NextcladeTask as nextclade_task {
        input:
        ha_fasta = irma_task.ha_fasta,
        na_fasta = irma_task.na_fasta
    }

    call tasks.NextcladeParse as nextclade_parese_task {
        input:
            HA_tsv=nextclade_task.HA_nextclade,
            NA_tsv=nextclade_task.NA_nextclade,
            report=blast_parse_task.report_txt,
            sample_name=sample_name
    }

    output {
        File report=blast_parse_task.report_txt
    }
}






