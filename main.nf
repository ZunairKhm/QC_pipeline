    #!/usr/bin/env nextflow
    nextflow.enable.dsl=2

    include { FASTP } from './modules/nf-core/fastp/main'
    include { SYLPH_PROFILE } from './modules/nf-core/sylph/profile/main' 
    include { SYLPH_QUERY } from './modules/nf-core/sylph/query/main'               

    params.kingfisher   = null
    params.filepath       = null
    params.outdir         = 'results'


    workflow {

        reads_ch = Channel.empty()

        if (params.kingfisher && !params.filepath) {

            def acc_file = file(params.kingfisher)

            accessions_ch = acc_file.exists()
                ? Channel.fromPath(acc_file)
                        .splitCsv(header: false)
                        .map { row -> tuple([id: row[0]], row[0]) }
                : Channel.from(params.kingfisher.split(','))
                        .map { acc -> tuple([id: acc], acc) }

            KINGFISHER_GET(accessions_ch)

            reads_ch = KINGFISHER_GET.out.reads
                .map { meta, fq_list ->
                    fq_list = fq_list instanceof List ? fq_list : [fq_list]
                    fq_list.sort()
                    tuple([id: meta.id, single_end: fq_list.size() == 1], fq_list)
                }

    } else if (params.filepath && !params.kingfisher) {
        

            def fastq_dir = file(params.filepath)

            reads_ch = Channel
                .fromPath(fastq_dir + '/*.fastq.gz')
                .ifEmpty { error "No FASTQ files found in directory: ${fastq_dir}" }
                .map { f ->
                    def matcher = f.name =~ /(.+?)(?:[_-]R?([12])(?:_001)?)?\.fastq\.gz/
                    if (!matcher) error "File '${f.name}' does not match FASTQ pattern"
                    def id = matcher[0][1]
                    def is_single_end = matcher[0][2] == null
                    tuple([id: id, single_end: is_single_end], f)
                }
                .groupTuple(by: 0)
                .map { meta, files -> tuple(meta, files.sort()) }

           

        }
    else {
        error "No valid input source specified. Use either --kingfisher and a file of line separated accessions or list accessions directly, or --filepath and a file path to a directory containing fastq files."
    }

    reads_ch.view { meta, files -> "Raw file recieved: ${meta.id} -> ${files*.name}" }

    if (!params.deacon_index) {
        error "Please provide --deacon_index"
    }

    deacon_index_ch = Channel.value(file(params.deacon_index))

    // Run DEACON
    DEACON_FILTER(reads_ch, deacon_index_ch)
    //deacon_ch.filtered_reads.view { meta, files -> "Deacon filtered: ${meta.id} -> ${files*.name}" }

    // Run FASTP
        FASTP(
            DEACON_FILTER.out.filtered_reads,
            [],  // adapter_fasta
            [],  // discard_trimmed_pass
            [],  // save_trimmed_fail
            [] // save_merged
        )

    //fastp_ch.reads.view { meta, files -> "FASTP output: ${meta.id} -> ${files*.name}" }

// Path to Sylph database
sylph_95_db_ch = Channel.value(file(params.sylph_95_db))
sylph_99_db_ch = Channel.value(file(params.sylph_99_db))

// Run Sylph
sylph_profile_ch = SYLPH_PROFILE(FASTP.out.reads, sylph_95_db_ch)
sylph_query_ch = SYLPH_QUERY(FASTP.out.reads, sylph_99_db_ch)

    // View outputs
    sylph_profile_ch.profile_out.view { meta, tsv_files ->
    "Sylph profile: ${meta.id} -> ${tsv_files*.name}"
    }
    sylph_query_ch.query_out.view { meta, tsv_files ->
    "Sylph query: ${meta.id} -> ${tsv_files*.name}"
    }

    } // end workflow
       
    


    // ----------- Processes --------------

    process KINGFISHER_GET {
        tag "$meta.id"
        //publishDir "${params.outdir}/raw_data", mode: 'symlink'

        input:
        tuple val(meta), val(accession)

        output:
        tuple val(meta), path("*.fastq.gz"), emit: reads    

        errorStrategy 'retry'
        maxRetries 10
        maxErrors '-1' 

        script:
        """
        kingfisher get -r ${accession} -m ena-ascp ena-ftp aws-http --output-directory . -t $task.cpus -f "fastq.gz"
        """
    }

process DEACON_FILTER {
    tag "$meta.id"
    publishDir "${params.outdir}/deacon_log",  mode: 'copy', pattern: "*.deacon.log"
    //publishDir "${params.outdir}/processed_reads", mode: 'copy', pattern: "*filt*.fq.gz" 

    input:
    tuple val(meta), path(fq_files)
    path(deacon_index_file)

    output:
    tuple val(meta), path("${meta.id}.filt*.fq.gz"), emit: filtered_reads  // More specific pattern
    tuple val(meta), path("${meta.id}.deacon.log"), emit: deacon_log

    script:
    def cmd
    def outputs

    if (fq_files.size() == 2) {
        // paired-end
        cmd = "deacon filter -d ${deacon_index_file} ${fq_files[0]} ${fq_files[1]} -o ${meta.id}.filt.R1.fq.gz -O ${meta.id}.filt.R2.fq.gz -s ${meta.id}.deacon.log"
        outputs = "${meta.id}.filt.R1.fq.gz ${meta.id}.filt.R2.fq.gz"
    } else {
        // single-end
        cmd = "deacon filter -d ${deacon_index_file} ${fq_files[0]} -o ${meta.id}.filt.fq.gz -s ${meta.id}.deacon.log"
        outputs = "${meta.id}.filt.fq.gz"
    }

    """
    echo "Processing sample: ${meta.id}"
    echo "Input files: ${fq_files*.name.join(' ')}"
    echo "Index file: ${deacon_index_file}"

    $cmd
    ls -lh $outputs
    """
}






