
workflow {
    ref = Channel.fromPath(params.ref)
    | collect
   
raw_alignments = Channel
    .fromPath(
        "/varidata/research/projects/laird/Ava/chrons_human_data/results/alignments/*/*.cram"
    )
    .map { cram ->
        def id = cram.baseName
        def crai = file("${cram}.crai")

        tuple(id, cram, crai)
    }
    .combine(ref)
    .map { row ->

        def id = row[0]
        def cram = row[1]
        def crai = row[2]

        def ref_files = row[3..-1]

        tuple(id, cram, crai, ref_files)
    }


raw_alignments | (trace_reads)

}



process trimmed_fastqc {
    container 'community.wave.seqera.io/library/biscuit_dupsifter_fastqc_samtools:3899cf51dd7f8705'
    cpus 2
    publishDir "./results/fastqc/trimmed", mode: 'copy', overwrite: true
    shell = ["/bin/bash", "-euo", "pipefail"]
    errorStrategy = 'ignore'

    input:
        tuple val(id), path(cram), path(crai), path(ref)

    output:
        tuple val(id), path("*.zip"), path("*.html")

    script:
    """
        mkfifo read_1 read_2

        cat <(samtools view -hT ${ref[1]} ${cram} | awk '\$1 ~ /^@/ || (\$2 != 165 && \$2 != 101)') \
            <(samtools view -T ${ref[1]} --require-flags PAIRED,UNMAP,MREVERSE,READ2 \
                --exclude-flags PROPER_PAIR,MUNMAP,REVERSE,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY,READ1 \
                --add-flags REVERSE ${cram}) \
            <(samtools view -T ${ref[1]} --require-flags PAIRED,UNMAP,MREVERSE,READ1 \
                --exclude-flags PROPER_PAIR,MUNMAP,REVERSE,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY,READ2 \
                --add-flags REVERSE ${cram}) \
            | samtools fastq -1 read_1 -2 read_2 & pid1=\$!
        
        fastqc --nogroup stdin:${id}_R1 < read_1 & pid2=\$!
        fastqc --nogroup stdin:${id}_R2 < read_2 & pid3=\$!

        wait \$pid1 \$pid2 \$pid3
        rm read_1 read_2
    """
}

process trace_reads {

    container 'community.wave.seqera.io/library/biscuit_dupsifter_fastqc_samtools:3899cf51dd7f8705'
    cpus 2

    publishDir "./results/read_tracing", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(cram), path(crai), path(ref)

    output:
        tuple val(id),
            path("${id}.deduplicated.bam"),
            path("${id}.deduplicated.bam.bai"),
            path("${id}_stats.tsv")

    script:
    """
        mapping_stats() {
            local bam_file="\$1"
            local view_flags="\$2"
            local label="\$3"

            local reads
            reads=\$(samtools view -c \$view_flags "\$bam_file")

            local bases
            bases=\$(samtools view \$view_flags "\$bam_file" |
                awk '\$10 != "*" {s+=length(\$10)} END {print s}')

            local aln_bases
            aln_bases=\$(samtools view -h \$view_flags "\$bam_file" |
                samtools depth -g 0xFFF - |
                awk '{s+=\$3} END {print s}')

            local cov_bases
            cov_bases=\$(samtools view -h \$view_flags "\$bam_file" |
                samtools coverage -H - |
                awk '{s+=\$5} END {print s}')

            printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \
                "\$label" "\$reads" "\$bases" "\$aln_bases" "\$cov_bases"
        }


        # ============================================================
        # CRAM → COLLATE → DUPSIFTER → SORT
        # ============================================================

        samtools view -hbT ${ref[1]} ${cram} \
            | samtools collate \
                -T collate_${id} \
                --output-fmt sam \
                -O - \
            | dupsifter ${ref[1]} \
            | samtools sort \
                -o ${id}.deduplicated.bam


        # ============================================================
        # BISCUIT FILTERING
        # ============================================================

        biscuit bsconv -y 0.01 \
            ${ref[1]} \
            ${id}.deduplicated.bam \
            ${id}.bsconv.bam

        mv ${id}.bsconv.bam ${id}.deduplicated.bam


        # ============================================================
        # INDEX FINAL BAM
        # ============================================================

        samtools index ${id}.deduplicated.bam


        # ============================================================
        # STATS FROM FINAL BAM
        # ============================================================

        cat <<EOF > ${id}_stats.tsv
label\\tn_reads\\tn_bases_in_bam\\tn_ref_aln_bases\\tn_bases_covered
\$(mapping_stats "${id}.deduplicated.bam" "" "Bisulfite_Conversion_Filter")
\$(mapping_stats "${id}.deduplicated.bam" "-F 0x900" "Primary")
\$(mapping_stats "${id}.deduplicated.bam" "-F 0x904 -q 41" "MAPQ > 40")
\$(mapping_stats "${id}.deduplicated.bam" "-F 0xd04 -q 41" "Remove_Duplicates")
        EOF
    """
}