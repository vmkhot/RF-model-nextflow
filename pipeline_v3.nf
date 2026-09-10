#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.out_dir         = "output_dir_1"
params.clusters_file   = null
params.fasta_file      = null
params.biomes_file     = null
params.taxa_file       = null
params.num_seqs        = 25
params.min_seqs_own_model = 30
params.wanted_taxonomy = null
params.biome1          = null
params.biome2          = null
params.help            = false
params.do_trimal       = false

if (params.help) {
    log.info """
    Usage:
        nextflow run pipeline.nf --clusters_file file.tsv --fasta_file file.faa --biomes_file file.tsv --taxa_file file.tsv --wanted_taxonomy TAX --biome1 "string" --biome2 "string"

    Options:
        --out_dir           Output directory (default: output_dir_1)
        --clusters_file     Clusters .tsv file
        --fasta_file        Input protein FASTA file
        --biomes_file       Biomes file
        --taxa_file         Taxonomy file
        --num_seqs          Number of sequences to sample per cluster (default: 25)
        --min_seqs_own_model Minimum number of sequences in a cluster for it to form its own model (otherwise go to mixed)
        --wanted_taxonomy   Taxonomy filter
        --biome1            Biome filter 1
        --biome2            Biome filter 2
    """
    exit 0
}

// Process: clusters_to_fasta
process CLUSTERS_TO_FASTA {
    cpus 5
    memory '100 GB'
    time '3h'

    publishDir "${params.out_dir}/cluster_fastas", mode: 'copy'

    output:
    path "*.faa", emit: faa
    path "clusters.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

    python3 /vast/ri65fin/mcp_structures/nextflow/bin/clusters_to_fasta.py \\
        -clusters_file ${params.clusters_file} \\
        -fasta_file ${params.fasta_file} \\
        -biomes_file ${params.biomes_file} \\
        -taxa_file ${params.taxa_file} \\
        --num_seqs ${params.num_seqs} \\
        --min_seqs_own_model ${params.min_seqs_own_model} \\
        --out_dir \${PWD} \\
        --wanted_taxonomy "${params.wanted_taxonomy}" \\
        --biome1 "${params.biome1}" \\
        --biome2 "${params.biome2}" > clusters.log 2>&1
    """
}

// Process: foldseek
process FOLDSEEK {
    cpus 30
    memory '50 GB'
    time '1h'
    maxForks 4      // Only open 4 processes

    errorStrategy 'retry'   // retry if exit code =! 0
    maxRetries 2            // retry max 2x

    publishDir "${params.out_dir}/foldseek", mode: 'copy'

    queue 'gpu-veo,gpu-test,gpu'
    clusterOptions '--gres=gpu:1'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_db", emit: db
    path "foldseek.log", emit: log


    script:
    """
    module load nvidia/cuda/12.6.0
    export PATH=/home/groups/VEO/tools/foldseek/v10-941cd33/bin/:\$PATH

    db_dir=${fasta_file.baseName}_db
    mkdir -p \$db_dir

    foldseek createdb ${fasta_file} \$db_dir/${fasta_file.baseName} \\
        --prostt5-model /vast/ri65fin/mcp_structures/MCP_struct/prostt5_10_941_weights \\
        --gpu 1 > foldseek.log 2>&1
    """
}

// Process: foldmason
process FOLDMASON {
    cpus 30
    memory '100 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/foldmason", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path db_folder

    output:
    path "${db_folder}_out", emit: db
    path "foldmason.log", emit: log


    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate foldmason

    # create output folder for this cluster
    mkdir -p ${db_folder}_out

    # extract basename inside bash
    cluster_name=\$(basename ${db_folder} _db)

    # run Foldmason
    foldmason structuremsa ${db_folder}/\$cluster_name ${db_folder}_out/\${cluster_name}_MSTA > foldmason.log 2>&1 
    """
}

process TRIMAL {
    cpus 4
    memory '10 GB'
    time '1h'

    publishDir "${params.out_dir}/RF_results/foldmason", mode: 'copy'

    input:
    path foldmason_out

    output:
    path "${foldmason_out.baseName}_trimmed", emit: faa
    path "trimal.log", emit: log

    script:
    """
    clust_num=\$(basename ${foldmason_out} _db_out)

    outdir=${foldmason_out.baseName}_trimmed
    mkdir -p \$outdir

    /home/groups/VEO/tools/trimal/v1.5.0/trimal/source/trimal -in ${foldmason_out}/\${clust_num}_MSTA_aa.fa \\
           -out \$outdir/\${clust_num}_MSTA_trimmed.fa \\
           -gappyout \\
           -colnumbering > trimal.log 2>&1
    """
}

process RF_MODEL {
    cpus 30
    memory '100 GB'
    time '6h'

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/RF_results/raw", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path msa_folder

    output:
    path "RF_results/*_Feature_importance.tsv" , emit: features
    path "RF_results/*_RandomForest_Report.tsv", emit: reports
    path "RF_results/*_data_df.tsv", emit: data_df
    path "rf_model.log", emit: log


    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

    clust_num=\$(basename ${msa_folder} | sed 's/_db_out.*//')

    # decide input file: trimmed or untrimmed
    if [ -f ${msa_folder}/\${clust_num}_MSTA_trimmed.fa ]; then
        msa_file=${msa_folder}/\${clust_num}_MSTA_trimmed.fa
    else
        msa_file=${msa_folder}/\${clust_num}_MSTA_aa.fa
    fi

    python3 /vast/ri65fin/mcp_structures/nextflow/bin/rf_model.py \\
        "\$msa_file" \\
        ${params.biomes_file} >> rf_model.log 2>&1
    """
}

process MERGE_RF_RESULTS {
    publishDir "${params.out_dir}/RF_results", mode: 'copy'

    input:
    path feature_files
    path report_files

    output:
    path "Feature_importance.tsv"
    path "RandomForest_Report.tsv"

    script:
    """
    # Merge feature importance files
    awk 'FNR==1 && NR!=1 { next } { print }' \$(ls *_Feature_importance.tsv | sort) > Feature_importance.tsv

    # Merge report files
    awk 'FNR==1 && NR!=1 { next } { print }' \$(ls *_RandomForest_Report.tsv | sort) > RandomForest_Report.tsv

    """
}

// Workflow
workflow {
    // Call CLUSTERS_TO_FASTA, take only faa output, flatten
    fasta_ch = CLUSTERS_TO_FASTA().faa | flatten

    // Call FOLDSEEK with fasta_ch, take db output, flatten
    db_ch    = FOLDSEEK(fasta_ch).db | flatten

    // Call FOLDMASON with db_ch, take db output, flatten
    mason_ch = FOLDMASON(db_ch).db | flatten

    // Optionally call TRIMAL, otherwise pass through mason_ch
    msa_ch   = params.do_trimal ? TRIMAL(mason_ch).faa : mason_ch

    // Call RF_MODEL with msa_ch
    rf_ch = RF_MODEL(msa_ch)

    // Merge RF results
    MERGE_RF_RESULTS(
        rf_ch.features.collect(),
        rf_ch.reports.collect()
    )
}

