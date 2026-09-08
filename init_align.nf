
workflow {

    fastqs = Channel.fromFilePairs(
        '/varidata/research/projects/laird/Ava/chrons_human_data/raw_data/*_R{1,2}.fastq.gz',
        flat: true
    )

    ref = Channel.fromPath(params.ref)
    | collect  
    | map { refs -> tuple(refs) }

    fastqs
        .combine(ref)
        .map { row ->

        def id = row[0]

       def fastqs = row[1..2]   // ALWAYS explicit list

        def ref_fa = row[3..-1]            // ONLY the FASTA

        tuple(id, fastqs, ref_fa)
    }
   | trim_align
}
process trim_align {

    container 'community.wave.seqera.io/library/biscuit_cutadapt_samtools:15408f8efe67445f'
    cpus 64
    memory 150.GB

    publishDir "./results/alignments/${sample_id}",
        mode: 'copy',
        pattern: "*.cram*"

    publishDir "./results/trimming_reports/${sample_id}",
        mode: 'copy',
        pattern: "*.json"

    shell = ["/bin/bash", "-euo", "pipefail"]

    input:
        tuple val(sample_id), path(fastqs), path(ref)

    output:
        tuple val(sample_id), path("*.cram"), path("*.crai"), emit: alignments
        tuple val(sample_id), path("*.json"), emit: cutadapt_reports

    script:
        def id = sample_id

        """
        cutadapt --cores 8 --interleaved \
            --json ${id}_cutadapt_report.json \
            --nextseq-trim 20 --overlap 1 -n 7 \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC ${fastqs} \
            | paste - - - - - - - - \
            | awk -v 'FS=\\t' -v 'OFS=\\n' '{print \$5,\$6,\$7,\$8,\$1,\$2,\$3,\$4}' \
            | biscuit align -@ 60 -b 1 -p \
                -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' \
                ${ref[1]} - \
            | samtools sort \
                --write-index \
                --reference ${ref[1]} \
                -o ${id}.cram -
        """
}