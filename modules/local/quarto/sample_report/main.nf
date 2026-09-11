process QUARTO_SAMPLE_REPORT {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    // No container declared: the only published quarto+R image
    // (community.wave.seqera.io/library/quarto_..._pruned:e6e33e4fe2c03435) is a
    // *pruned* Wave build whose dart-sass can't compile the HTML theme
    // ("Theme file compilation failed"). This module therefore runs via conda -
    // use the conda_singularity profile.

    input:
    tuple val(meta), path(bam), path(bai), path(resistance), path(formatted), path(assessed) //sorted bam (empty for fasta input), resistance + formatted from MERGE_FORMAT, assessed.tsv (fasta input only, empty otherwise)
    path(mut_bed) //this is the bed file listing the locations of the mutations for the report plot
    val(protocol) //one of amplicon or wgs
    val(input_type) //one of short-read, long-read or fasta
    val(report_author) //name that you want to appear at the top of the report
    val(bibliography) //.bib .bibtex format
    val(database_name) //name of the database used
    val(cmv_resis_db_version) //date
    val(min_depth) //sequencing depth threshold for calling mutations
    val(cmv_merlin_name) //reference sequence name to list on report
    val(gene_positions) //this is the start/end positions of the AMR genes for plotting coverage.

    output:
    tuple val(meta), path("${meta.id}.html")
    tuple val(meta), path("${meta.id}.pdf"), optional: true, emit: pdf
    tuple val("${task.process}"), val('quarto'), eval('quarto -v'), emit: versions_quarto, topic: versions

    script:
    def ul27_pos = gene_positions['ul27']
    def ul51_pos = gene_positions['ul51']
    def ul54_pos = gene_positions['ul54']
    def ul56_pos = gene_positions['ul56']
    def ul89_a_pos = gene_positions['ul89_a']
    def ul89_b_pos = gene_positions['ul89_b']
    def ul97_pos = gene_positions['ul97']

    """
    cp ${projectDir}/bin/sample_report.qmd .
    cp ${bibliography} references.bib

cat > report-params.yml <<END_PARAMS
report_author: "${report_author}"
database_name: "${database_name}"
database_version: "${cmv_resis_db_version}"
sample: "${meta.id}"
single_end: "${meta.single_end}"
resistance: "${resistance}"
formatted: "${formatted}"
assessed: "${assessed}"
bam: "${bam}"
mut_bed: "${mut_bed}"
protocol: "${protocol}"
input_type: "${input_type}"
threshold: "${min_depth}"
reference: "${cmv_merlin_name}"
ul27_pos: "${ul27_pos}"
ul51_pos: "${ul51_pos}"
ul54_pos: "${ul54_pos}"
ul56_pos: "${ul56_pos}"
ul89_a_pos: "${ul89_a_pos}"
ul89_b_pos: "${ul89_b_pos}"
ul97_pos: "${ul97_pos}"
END_PARAMS

    quarto render sample_report.qmd \\
        --to html \\
        --execute-params report-params.yml \\
        --output ${meta.id}.html

    # PDF (typst) is best-effort: a failed render must not fail the pipeline.
    if quarto render sample_report.qmd --to typst --execute-params report-params.yml --output ${meta.id}.pdf ; then
        echo "PDF report written: ${meta.id}.pdf"
    else
        echo "WARNING: PDF (typst) render failed for ${meta.id} - the HTML report was still produced." >&2
    fi
    """
}

