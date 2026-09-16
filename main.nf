    #!/usr/bin/env nextflow
    nextflow.enable.dsl=2

    include { FASTP } from './modules/nf-core/fastp/main'
    include { SYLPH_PROFILE } from './modules/nf-core/sylph/profile/main' 
    include { SYLPH_QUERY } from './modules/nf-core/sylph/query/main'     
    include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'          

    params.kingfisher   = null
    params.filepath       = null
    params.outdir         = 'results'
    params.runkraken2     = false
    params.runsylphquery  = false
    params.runtracsalign  = false


    params.filter_method = 'nohuman' // Options: 'deacon' or 'nohuman' or 'none'

    // Set --run_label when running batches into a shared --outdir, so each batch writes its
    // own manifest instead of overwriting (and wiping) the previous batch's.
    params.run_label = null

    def runSuffix() { params.run_label ? ".${params.run_label}" : '' }

    // Each stage drops a small file of the sample IDs that reached it into this directory;
    // the onComplete handler reads them back to build the completion manifest.
    def manifestDir() { "${params.outdir}/.manifest${runSuffix()}" }

    def sampleIdsFrom(String name) {
        def f = file("${manifestDir()}/${name}")
        f.exists() ? f.readLines()*.trim().findAll { it } as Set : [] as Set
    }

    // Failed samples are skipped (errorStrategy 'ignore'), so compare the samples that
    // entered the pipeline against those each stage actually produced output for.
    // Must live in a function: `params`/`log` don't resolve inside the onComplete closure.
    def writeSampleManifest() {
        def all = sampleIdsFrom('all_samples.txt')
        if (!all) {
            log.warn "No samples were recorded - skipping completion manifest."
            return
        }

        def stages = [
            ['fastp', sampleIdsFrom('fastp_samples.txt')],
            ['sylph_profile', sampleIdsFrom('sylph_profile_samples.txt')]
        ]
        if (params.runsylphquery) stages << ['sylph_query', sampleIdsFrom('sylph_query_samples.txt')]
        if (params.runkraken2)    stages << ['kraken2', sampleIdsFrom('kraken2_samples.txt')]
        if (params.runtracsalign) stages << ['tracs_align', sampleIdsFrom('tracs_align_samples.txt')]

        def report = file("${params.outdir}/sample_manifest${runSuffix()}.tsv")
        report.parent.mkdirs()
        report.text = (
            [(['sample_id'] + stages.collect { stage -> stage[0] }).join('\t')] +
            all.sort().collect { id ->
                ([id] + stages.collect { stage -> stage[1].contains(id) ? 'PASS' : 'FAIL' }).join('\t')
            }
        ).join('\n') + '\n'

        log.info "Sample manifest: ${report}"
        stages.each { stage ->
            def missing = all - stage[1]
            if (missing) {
                log.warn "${missing.size()}/${all.size()} sample(s) produced no ${stage[0]} output: ${missing.sort().join(', ')}"
            }
        }
    }

    workflow {

        // Clear last run's stage files: collectFile writes nothing for a stage that produced
        // no output, so stale files would otherwise be read back as false PASSes.
        file(manifestDir()).deleteDir()

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
    reads_ch.map { meta, _files -> meta.id }
        .collectFile(name: 'all_samples.txt', newLine: true, storeDir: manifestDir())

// ---------------------------------------------------------
    // Filtering Logic (Deacon vs Nohuman)
    // ---------------------------------------------------------
    def filtered_reads_ch = Channel.empty()

    if (params.filter_method == 'deacon') {
        if (!params.deacon_index) {
            error "Please provide --deacon_index when using filter_method 'deacon'"
        }
        deacon_index_ch = Channel.value(file(params.deacon_index))
        
        DEACON_FILTER(reads_ch, deacon_index_ch)
        filtered_reads_ch = DEACON_FILTER.out.filtered_reads

    } else if (params.filter_method == 'nohuman') {
        if (!params.nohuman_db) {
            error "Please provide --nohuman_db when using filter_method 'nohuman'"
        }
        nohuman_db_ch = Channel.value(file(params.nohuman_db))
        
        NOHUMAN_FILTER(reads_ch, nohuman_db_ch)
        filtered_reads_ch = NOHUMAN_FILTER.out.filtered_reads

    } else if (params.filter_method == 'none') {
        filtered_reads_ch = reads_ch
    } else {
        error "Invalid --filter_method specified. Must be 'deacon' or 'nohuman'."
    }

    // Run FASTP
        FASTP(
            filtered_reads_ch,
            [],  // adapter_fasta
            [],  // discard_trimmed_pass
            [],  // save_trimmed_fail
            [] // save_merged
        )

    //fastp_ch.reads.view { meta, files -> "FASTP output: ${meta.id} -> ${files*.name}" }
    FASTP.out.reads.map { meta, _reads -> meta.id }
        .collectFile(name: 'fastp_samples.txt', newLine: true, storeDir: manifestDir())

// Path to Sylph database
sylph_95_db_ch = channel.value(file(params.sylph_95_db))
sylph_99_db_ch = channel.value(file(params.sylph_99_db))

// Run Sylph
sylph_profile_ch = SYLPH_PROFILE(FASTP.out.reads, sylph_95_db_ch)
if (params.runsylphquery) {    
sylph_query_ch = SYLPH_QUERY(FASTP.out.reads, sylph_99_db_ch)
}

    // View outputs
    sylph_profile_ch.profile_out.view { meta, tsv_files ->
    "Sylph profile: ${meta.id} -> ${tsv_files*.name}"
    }
    sylph_profile_ch.profile_out.map { meta, _tsv_files -> meta.id }
        .collectFile(name: 'sylph_profile_samples.txt', newLine: true, storeDir: manifestDir())

    if (params.runsylphquery) {
        sylph_query_ch.query_out.view { meta, tsv_files ->
        "Sylph query: ${meta.id} -> ${tsv_files*.name}"
        }
        sylph_query_ch.query_out.map { meta, _tsv_files -> meta.id }
            .collectFile(name: 'sylph_query_samples.txt', newLine: true, storeDir: manifestDir())
    }


if (params.runkraken2) {
    KRAKEN2_KRAKEN2 (
            FASTP.out.reads,
            file(params.kraken2_db),
            false,                                           // Don't save classified FASTQ files
            false                                            // Don't save the massive raw assignment file
        )

    KRAKEN2_KRAKEN2.out.report.map { meta, _report -> meta.id }
        .collectFile(name: 'kraken2_samples.txt', newLine: true, storeDir: manifestDir())
}

if (params.runtracsalign) {
    if (!params.tracs_refseqs) {
        error "Please provide --tracs_refseqs (path to reference fasta file(s)) when using --runtracsalign"
    }
    if (!params.tracs_database) {
        error "Please provide --tracs_database (a sourmash .sbt.zip) when using --runtracsalign"
    }
    tracs_refseqs_ch  = Channel.value(file(params.tracs_refseqs))
    tracs_database_ch = Channel.value(file(params.tracs_database))

    TRACS_ALIGN(FASTP.out.reads, tracs_refseqs_ch, tracs_database_ch)

    TRACS_ALIGN.out.align_out.view { meta, dir -> "TRACS align: ${meta.id} -> ${dir.name}" }
    TRACS_ALIGN.out.align_out.map { meta, _dir -> meta.id }
        .collectFile(name: 'tracs_align_samples.txt', newLine: true, storeDir: manifestDir())
}


    workflow.onComplete {
        writeSampleManifest()
    }

    } // end workflow




    // ----------- Processes --------------

    process KINGFISHER_GET {
        tag "$meta.id"
        //publishDir "${params.outdir}/raw_data", mode: 'copy'

        // Retry transient download failures, then skip an accession that keeps failing
        // (withdrawn/restricted). A plain 'retry' aborts the entire run once retries run out.
        errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
        maxRetries 3

        input:
        tuple val(meta), val(accession)

        output:
        tuple val(meta), path("*.fastq.gz"), emit: reads

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

process NOHUMAN_FILTER {
    tag "$meta.id"
    publishDir "${params.outdir}/nohuman_log", mode: 'copy', pattern: "*.report"

    input:
    tuple val(meta), path(fq_files)
    path(nohuman_db)

    output:
    // nohuman natively appends 'nohuman' to the filename before .fastq.gz
    tuple val(meta), path("*nohuman*.fq.gz"), emit: filtered_reads
    tuple val(meta), path("*.report"), emit: nohuman_report

    script:
    def cmd
    if (fq_files.size() == 2) {
        cmd = "nohuman -D ${nohuman_db} -k ${meta.id}.kraken.out -r ${meta.id}.kraken.report ${fq_files[0]} ${fq_files[1]}"
    } else {
        cmd = "nohuman -D ${nohuman_db} -k ${meta.id}.kraken.out -r ${meta.id}.kraken.report ${fq_files[0]}"
    }

    """
    echo "Processing sample: ${meta.id} with nohuman"
    echo "Database: ${nohuman_db}"

    $cmd
    """
}

process TRACS_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/tracs_align", mode: 'copy'

    input:
    tuple val(meta), path(fq_files)
    path(tracs_refseqs)
    path(tracs_database)

    output:
    tuple val(meta), path("${meta.id}"), emit: align_out

    script:
    def input_reads = fq_files instanceof List ? fq_files.join(' ') : fq_files
    // $HOME is read-only in the container, so point matplotlib/fontconfig at the work dir
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"

    tracs align \\
        -i ${input_reads} \\
        --database ${tracs_database} \\
        --refseqs ${tracs_refseqs} \\
        -o ${meta.id} \\
        -p ${meta.id} \\
        -t $task.cpus
    """
}


