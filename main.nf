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

    // reads a "<unit>\t<v1>\t<v2>..." file from the manifest dir into unit -> column value
    def unitColumnFrom(String name, int col) {
        def f = file("${manifestDir()}/${name}")
        def out = [:]
        if (f.exists()) {
            f.readLines().each { line ->
                def parts = line.trim().split('\t')
                if (parts.size() > col) out[parts[0]] = parts[col]
            }
        }
        return out
    }

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

        // units are merged from one or more runs, so report that alongside the stage outcomes
        def unitRuns   = unitColumnFrom('unit_runs.txt', 1)
        def mergeState = unitColumnFrom('merge_status.txt', 3)

        def report = file("${params.outdir}/sample_manifest${runSuffix()}.tsv")
        report.parent.mkdirs()
        report.text = (
            [(['sample_id', 'runs', 'merge'] + stages.collect { stage -> stage[0] }).join('\t')] +
            all.sort().collect { id ->
                ([id, unitRuns.get(id, '1'), mergeState.get(id, 'single_run')] +
                 stages.collect { stage -> stage[1].contains(id) ? 'PASS' : 'FAIL' }).join('\t')
            }
        ).join('\n') + '\n'

        def partial = mergeState.findAll { _k, v -> v == 'PARTIAL' }.keySet()
        if (partial) {
            log.warn "${partial.size()} unit(s) merged from fewer runs than expected: ${partial.sort().join(', ')}"
        }

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

            // kingfisher annotate needs a file, so materialise a comma-separated list into one
            accession_file_ch = acc_file.exists()
                ? channel.fromPath(acc_file)
                : channel.of(params.kingfisher.split(',') as List).flatten()
                         .collectFile(name: 'accessions.txt', newLine: true)

            KINGFISHER_ANNOTATE(accession_file_ch)

            // One unit per sample. Runs of a sample are grouped by sample_accession; where a
            // sample mixes paired and single runs the single ones are dropped, since the two
            // layouts cannot be concatenated into one pair of files.
            sample_rows_ch = KINGFISHER_ANNOTATE.out.tsv
                .splitCsv(header: true, sep: '\t')
                .map { row -> tuple(row.sample_accession, row) }
                .groupTuple()

            units_ch = sample_rows_ch.flatMap { sample, rows ->
                def paired = rows.findAll { (it.library_layout ?: '').toUpperCase() != 'SINGLE' }
                def keep   = paired ?: rows
                def single_end = paired.isEmpty()
                // name the unit after its first run accession; sorted so the id is stable
                // across runs, otherwise -resume would miss every merged unit
                def id = keep.collect { it.run }.sort().first()
                keep.collect { r ->
                    tuple(groupKey([id: id, single_end: single_end, n_runs: keep.size()], keep.size()), r.run)
                }
            }

            units_ch
                .map { key, run -> "${key.id}\t${run}\t${key.single_end}\t${key.n_runs}" }
                .collectFile(name: "units${runSuffix()}.tsv", newLine: true, storeDir: "${params.outdir}/metadata",
                             seed: "unit_id\trun\tsingle_end\tn_runs")

            // record what the mixed-layout rule discarded, so the loss is visible
            sample_rows_ch
                .flatMap { sample, rows ->
                    def paired = rows.findAll { (it.library_layout ?: '').toUpperCase() != 'SINGLE' }
                    paired ? rows.findAll { (it.library_layout ?: '').toUpperCase() == 'SINGLE' }
                                 .collect { "${sample}\t${it.run}\tsingle_end_run_in_mixed_sample" } : []
                }
                .collectFile(name: "dropped_runs${runSuffix()}.tsv", newLine: true, storeDir: "${params.outdir}/metadata",
                             seed: "sample_accession\trun\treason")

            units_ch
                .map { key, _run -> "${key.id}\t${key.n_runs}" }
                .unique()
                .collectFile(name: 'unit_runs.txt', newLine: true, storeDir: manifestDir())

            KINGFISHER_GET(units_ch)

            // groupKey releases each unit as soon as its n_runs downloads land, rather than
            // waiting for the whole accession list. remainder keeps partly-failed units.
            grouped_ch = KINGFISHER_GET.out.reads
                .groupTuple(by: 0, remainder: true)
                .map { key, per_run -> tuple([id: key.id, single_end: key.single_end, n_runs: key.n_runs], per_run) }

            // Only a unit that was *meant* to have one run skips merging. A multi-run unit
            // reduced to one surviving download must still go through MERGE_RUNS, otherwise it
            // would pass through looking complete with no PARTIAL record.
            branched_ch = grouped_ch.branch { meta, per_run ->
                one_run: per_run.size() == 1 && meta.n_runs == 1
                many_runs: true
            }

            single_run_ch = branched_ch.one_run
                .map { meta, per_run ->
                    def fq = per_run.flatten().sort { a, b -> a.name <=> b.name }
                    // ENA sometimes ships a third, unpaired file (<run>.fastq.gz) next to
                    // <run>_1/_2; nohuman only takes one or two inputs, so without this the
                    // orphan alone went through as a single-end sample. MERGE_RUNS already
                    // drops such orphans for multi-run units.
                    def mates = fq.findAll { it.name ==~ /.*_[12]\.fastq\.gz/ }
                    if (mates.size() == 2) fq = mates
                    tuple([id: meta.id, single_end: fq.size() == 1], fq)
                }

            MERGE_RUNS(branched_ch.many_runs.map { meta, per_run ->
                tuple(meta, per_run.flatten().sort { a, b -> a.name <=> b.name })
            })

            MERGE_RUNS.out.log
                .splitCsv(header: true, sep: '\t')
                .map { row -> "${row.unit_id}\t${row.expected_runs}\t${row.merged_runs}\t${row.status}" }
                .collectFile(name: 'merge_status.txt', newLine: true, storeDir: manifestDir())

            // a one-file glob yields a bare Path, not a list; downstream does fq_files.size()
            merged_ch = MERGE_RUNS.out.merged.map { meta, fq ->
                tuple([id: meta.id, single_end: meta.single_end], fq instanceof List ? fq : [fq])
            }

            reads_ch = single_run_ch.mix(merged_ch)

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

    if (params.filter_method == 'none') {
        filtered_reads_ch = reads_ch
    } else if (params.filter_method == 'deacon') {
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

    } else  {
        error "Invalid --filter_method specified. Must be 'deacon', 'nohuman', or 'none'."
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

    // Reference selection is pulled out ahead of the alignment: it records how much work
    // each sample is (references x input GB) and its gather result is handed to the
    // alignment so the reads are not sketched and gathered a second time.
    TRACS_GATHER(FASTP.out.reads, tracs_database_ch)

    tracs_sized_ch = TRACS_GATHER.out.sized.map { meta, fq, csv, refcount, input_gb ->
        def rc = refcount as Integer
        tuple(meta, fq, csv, rc, rc * (input_gb as Double))
    }

    // the inputs to the runtime model, kept because work dirs get cleaned up long before
    // there is enough data to re-fit it
    TRACS_GATHER.out.sized
        .map { meta, _fq, _csv, refcount, input_gb -> "${meta.id}\t${refcount}\t${input_gb}" }
        .collectFile(name: "tracs_sizing${runSuffix()}.tsv", newLine: true,
                     storeDir: "${params.outdir}/metadata",
                     seed: "sample\trefcount\tinput_gb")

    TRACS_ALIGN(tracs_sized_ch, tracs_refseqs_ch, tracs_database_ch)

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
        tag "$accession"
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

// tracs align picks its references with sourmash and only then starts aligning. Running
// just the reference-selection half up front records how much work each sample is
// (references x input GB; before the read screen, runtime tracked that product at r=0.987
// across 759 aligns) and, since the patched tracs accepts a precomputed gather csv
// (TRACS_GATHER_CSV), the result is handed to TRACS_ALIGN rather than computed twice.
process TRACS_GATHER {
    tag "$meta.id"

    input:
    tuple val(meta), path(fq_files)
    path(tracs_database)

    output:
    // the gather csv goes to TRACS_ALIGN so it is not recomputed there
    tuple val(meta), path(fq_files), path("${meta.id}_sourmash_hits.csv"), env('REFCOUNT'), env('INPUT_GB'), emit: sized

    script:
    def fq_list     = fq_files instanceof List ? fq_files : [fq_files]
    def input_reads = fq_list.join(' ')
    def py_files    = fq_list.collect { "\"${it}\"" }.join(', ')
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TRACS_GATHER_CORES=\$(nproc)
    # fastgather (branchwater), as batches 00-06 used, rather than upstream's classic gather:
    # it can break ties differently and shift the reference set by one or two on diverse
    # samples, but keeps new batches consistent with the old ones and avoids classic gather's
    # 8-20 min (up to hours) single-threaded run. The image's tracs gathers classically
    # unless TRACS_FASTGATHER=1, and silently falls back to classic gather when the plugin
    # is missing, so its absence is made an error here instead.
    export TRACS_FASTGATHER=1
    python3 -c 'import sourmash_plugin_branchwater' || { echo "sourmash_plugin_branchwater missing - fastgather unavailable" >&2; exit 1; }

    # Calls tracs' own run_gather rather than reimplementing sketch+fastgather here. The
    # reference list is not the gather csv, it is that csv after the p_match / coverage-ratio
    # filter in utils.py - reimplementing the two commands alone counted 88 references on a
    # sample tracs aligned 28 of. Going through the library keeps this count identical to the
    # one align arrives at by construction, and keeps it that way if the filter changes.
    python3 - <<'PY'
import os, tempfile
from zipfile import ZipFile
from tracs.utils import run_gather, is_valid_sourmash_db

outdir = os.path.join(os.getcwd(), "")
temp_dir = os.path.join(tempfile.mkdtemp(dir=outdir), "")

# same database resolution align does before it gathers
db = "${tracs_database}"
if is_valid_sourmash_db(db):
    smdb = db
else:
    with ZipFile(db, "r") as archive:
        archive.extract("sourmashDB.sbt.zip", temp_dir)
        smdb = temp_dir + "sourmashDB.sbt.zip"

refs = run_gather(
    input_files=[${py_files}],
    databasefile=smdb,
    output=outdir + "${meta.id}_sourmash_hits",
    temp_dir=temp_dir,
)
open("refcount.txt", "w").write(str(len(refs)))
PY

    export REFCOUNT=\$(cat refcount.txt)
    # -L because the reads are symlinked in, so the link's own size is not the file's
    export INPUT_GB=\$(stat -Lc%s ${input_reads} | awk '{s+=\$1} END {printf "%.4f", s/1e9}')

    echo "sizing: \$REFCOUNT references, \$INPUT_GB GB input"
    """
}

process TRACS_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/tracs_align", mode: 'copy'

    input:
    tuple val(meta), path(fq_files), path(gather_csv), val(refcount), val(work_units)
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

    # Pawsey grants cpus in proportion to the memory request, not the literal `cpus`
    # directive - a 20GB request here is billed ~24 cores regardless of what's asked for.
    # \$task.cpus only ever bakes in the literal directive (verified empirically: it stays
    # at the requested value even when slurm inflates the real grant), and Nextflow's own
    # singularity wrapper strips SLURM_CPUS_PER_TASK from the container's environment
    # (it runs `env - PATH=... singularity exec ...`) - so neither reflects the true grant.
    # nproc reads the cgroup CPU affinity directly, which is unaffected by either problem.
    CORES=\$(nproc)

    # TRACS_GATHER already ran the gather (fastgather) on these reads, so its csv is handed
    # over and the sketch + gather are skipped here. TRACS_GATHER_CORES would only matter if
    # the csv were missing and align had to gather itself.
    export TRACS_GATHER_CORES=\$CORES
    export TRACS_GATHER_CSV="\$PWD/${gather_csv}"

    # One minimap2 pass against all ${refcount} matched references at once finds which
    # reads can align to which reference; each per-reference alignment then only sees
    # those reads instead of the whole read set (the alignment stage used to cost
    # references x reads). A read reported against more than MAX_HITS references is
    # aligned to all of them, so no reference ever misses a read it would have aligned.
    # The cap trades screen output volume against promiscuous reads: measured on the same
    # inputs, 200 beat 50 on a 39-reference sample (22 vs 25 min; 163 vs 143k promiscuous
    # pairs) while 50 beat 200 on a 261-reference sample (29 vs 40 min), so it scales as
    # ~8000/references, clamped to [50, 200]. Outputs are identical at any cap.
    export TRACS_SCREEN=1
    export TRACS_SCREEN_MAX_HITS=\$(( 8000 / ${refcount} > 200 ? 200 : (8000 / ${refcount} < 50 ? 50 : 8000 / ${refcount}) ))

    # Per-reference alignments are small once screened, so run many at once with few
    # threads each (samtools sort is capped at 256M/thread inside tracs, so 12 workers stay
    # well under 20 GB). The posterior stage (pure python, one reference at a time
    # upstream) runs in a process pool and is now the longest phase on reference-rich
    # samples (9 of 16.5 min on a 179-reference sample at 6 workers); each worker holds
    # ~1-2 GB for a large genome, hence the cap of 8 inside a 20 GB request.
    export TRACS_ALIGN_WORKERS=\$(( CORES / 2 > 2 ? CORES / 2 : 2 ))
    export TRACS_POSTERIOR_WORKERS=\$(( CORES < 8 ? CORES : 8 ))

    tracs align \\
        -i ${input_reads} \\
        --database ${tracs_database} \\
        --refseqs ${tracs_refseqs} \\
        -o ${meta.id} \\
        -p ${meta.id} \\
        -t \$CORES
    """
}

process KINGFISHER_ANNOTATE {
    tag "${accessions.name}"
    publishDir "${params.outdir}/metadata", mode: 'copy'

    // every unit is derived from this, so a silent skip here would empty the whole run
    errorStrategy 'terminate'

    input:
    path(accessions)

    output:
    path("run_metadata.tsv"), emit: tsv

    script:
    """
    kingfisher annotate --run-identifiers-list ${accessions} --all-columns -f tsv > run_metadata.tsv
    """
}

process MERGE_RUNS {
    tag "$meta.id"
    publishDir "${params.outdir}/merge_logs", mode: 'copy', pattern: "*.merge.tsv"

    input:
    tuple val(meta), path(fq_files)

    output:
    tuple val(meta), path("${meta.id}*.fq.gz"), emit: merged
    path "${meta.id}.merge.tsv"               , emit: log

    script:
    """
    # gzip members concatenate directly, so this needs no decompression
    shopt -s nullglob
    r1=( *_1.fastq.gz ); r2=( *_2.fastq.gz ); se=()
    for f in *.fastq.gz; do
        case "\$f" in *_1.fastq.gz|*_2.fastq.gz) ;; *) se+=( "\$f" );; esac
    done

    if [ "${meta.single_end}" = "true" ]; then
        if [ \${#r1[@]} -gt 0 ]; then
            echo "ERROR: ${meta.id} is single-end but paired files were downloaded" >&2; exit 1
        fi
        cat "\${se[@]}" > ${meta.id}.fq.gz
        merged=\${#se[@]}
    else
        if [ \${#r1[@]} -ne \${#r2[@]} ]; then
            echo "ERROR: ${meta.id} has \${#r1[@]} R1 but \${#r2[@]} R2 files" >&2; exit 1
        fi
        if [ \${#se[@]} -gt 0 ]; then
            echo "ERROR: ${meta.id} is paired-end but \${#se[@]} unpaired file(s) arrived" >&2; exit 1
        fi
        cat "\${r1[@]}" > ${meta.id}_1.fq.gz
        cat "\${r2[@]}" > ${meta.id}_2.fq.gz
        merged=\${#r1[@]}
    fi

    # a download skipped by errorStrategy 'ignore' would otherwise merge silently short
    status=COMPLETE
    [ "\$merged" -eq "${meta.n_runs}" ] || status=PARTIAL
    printf 'unit_id\\texpected_runs\\tmerged_runs\\tlayout\\tstatus\\n' > ${meta.id}.merge.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.n_runs}" "\$merged" \\
        "\$([ "${meta.single_end}" = "true" ] && echo single || echo paired)" "\$status" >> ${meta.id}.merge.tsv
    """
}
