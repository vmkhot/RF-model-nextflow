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

// // Process: clusters_to_fasta
// process CLUSTERS_TO_FASTA {
//     cpus 5
//     memory '100 GB'
//     time '3h'

//     publishDir "${params.out_dir}/cluster_fastas", mode: 'copy'

//     output:
//     path "*.faa", emit: faa
//     path "clusters.log", emit: log

//     script:
//     """
//     source /vast/groups/VEO/shared_data/Lydia_Varada/miniconda3/etc/profile.d/conda.sh
//     conda activate working-env

//     python3 /vast/ri65fin/mcp_structures/nextflow/bin/clusters_to_fasta.py \\
//         -clusters_file ${params.clusters_file} \\
//         -fasta_file ${params.fasta_file} \\
//         -biomes_file ${params.biomes_file} \\
//         -taxa_file ${params.taxa_file} \\
//         --num_seqs ${params.num_seqs} \\
//         --min_seqs_own_model ${params.min_seqs_own_model} \\
//         --out_dir \${PWD} \\
//         --wanted_taxonomy "${params.wanted_taxonomy}" \\
//         --biome1 "${params.biome1}" \\
//         --biome2 "${params.biome2}" > clusters.log 2>&1
//     """
// }

// Process: Mafft
process MAFFT {
    cpus 10
    memory '50 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/mafft", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path fasta_file

    output:
    path "${cluster_name}_mafft_aa.fa", emit: aligned
    path "mafft.log", emit: log

    script:
    cluster_name = fasta_file.baseName
    """
    mafft="/home/groups/VEO/tools/mafft/v7.505/bin/mafft"

    # create output folder for this cluster, matching TRIMAL/RF_MODEL naming
    mkdir -p ${cluster_name}_db_out

    # Fasta split in ref.fa (longest sequence) and remain.fa (what remains)
    awk -v reff="ref.fa" -v restf="remain.fa" '
    /^>/ {
        if (seq != "") {
            n++; headers[n] = header; seqs[n] = seq
        }
        header = \$0
        seq = ""
        next
    }
    {
        line = \$0
        gsub(/\\r/, "", line)
        gsub(/[ \\t]+\$/, "", line)
        seq = seq line
    }
    END {
        if (seq != "") { n++; headers[n] = header; seqs[n] = seq }
        maxlen = 0
        for (i = 1; i <= n; i++) {
            if (length(seqs[i]) > maxlen) { maxlen = length(seqs[i]); best = i }
        }
        print headers[best]"\\n"seqs[best] > reff
        for (i = 1; i <= n; i++) {
            if (i != best) {
                print headers[i]"\\n"seqs[i] >> restf
            }
        }
    }
    ' ${fasta_file}

    \$mafft --keeplength --add remain.fa ref.fa > ${cluster_name}_mafft_aa.fa 2> mafft.log
    """
}

process FASTTREE {
    cpus 10
    memory '50 GB'
    time '12h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/fasttree", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path msa

    output:
    path "${cluster_name}_combined_fastree_lg_pseudo.tree", emit: tree
    path "fasttree.log", emit: log

    script:
    cluster_name = msa.baseName.replaceAll(/_mafft_aa.*/, '')
    """
    fastTreeMP="/home/groups/VEO/tools/fastTreeMP/v2.1.11/FastTreeMP"

    \$fastTreeMP -lg -pseudo < ${msa} > ${cluster_name}_combined_fastree_lg_pseudo.tree 2> fasttree.log
    """
}

// process TRIMAL {
//     cpus 4
//     memory '10 GB'
//     time '1h'

//     publishDir "${params.out_dir}/RF_results/foldmason", mode: 'copy'

//     input:
//     path foldmason_out

//     output:
//     path "${foldmason_out.baseName}_trimmed", emit: faa
//     path "trimal.log", emit: log

//     script:
//     """
//     clust_num=\$(basename ${foldmason_out} _db_out)

//     outdir=${foldmason_out.baseName}_trimmed
//     mkdir -p \$outdir

//     /home/groups/VEO/tools/trimal/v1.5.0/trimal/source/trimal -in ${foldmason_out}/\${clust_num}_mafft_aa.fa \\
//            -out \$outdir/\${clust_num}_MSTA_trimmed.fa \\
//            -gappyout \\
//            -colnumbering > trimal.log 2>&1
//     """
// }

process RF_MODEL {
    cpus 30
    memory '100 GB'
    time '12h'

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/RF_results/raw", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path msa_file

    output:
    path "RF_results/*_Feature_importance.tsv" , emit: features
    path "RF_results/*_RandomForest_Report.tsv", emit: reports
    path "RF_results/*_data_df.tsv", emit: data_df
    path "rf_model.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

    python3 /vast/ri65fin/mcp_structures/nextflow/bin/rf_model_parallel.py \\
        "${msa_file}" \\
        ${params.biomes_file} >> rf_model.log 2>&1
    """
}

process PHYLOGLM {
    cpus 15
    memory '100 GB'
    time '12h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/phyloglm", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path data
    path feature
    path tree

    output:
    path "${cluster_name}_phyloglm_dataset.tsv", emit: phyloglm_dataset
    path "${cluster_name}_phyloglm.csv", emit: csv
    path "${cluster_name}_phyloglm_plot.png", emit: phyloglm_plot
    path "${cluster_name}_zvalue_plot.png", emit: zvalue_plot
    path "${cluster_name}_gini.tsv", emit: gini_dataset
    path "phyloglm.log", emit: log

    script:
    cluster_name = tree.baseName.replaceAll(/_combined_fastree_lg_pseudo.*/, '')
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate phyloglm

    # create_dataset_phyloglm.py: biome (data_df) + feature (feat_imp) -> dataset
    python3 /vast/ri65fin/mcp_structures/nextflow/bin/create_dataset_phyloglm.py \\
        ${feature} ${data} -o ${cluster_name}_phyloglm_dataset.tsv >> phyloglm.log 2>&1

    # phyloglm.R: dataset + tree -> csv results + plot
    Rscript /vast/ri65fin/mcp_structures/nextflow/bin/phyloglm.R \\
        ${cluster_name}_phyloglm_dataset.tsv ${tree} ${cluster_name}_gini.tsv\\
        ${cluster_name}_phyloglm.csv ${cluster_name}_phyloglm_plot.png >> phyloglm.log 2>&1

    # plot_zvalue_featimps.R: results + feature -> zvalue plot
    Rscript /vast/ri65fin/mcp_structures/nextflow/bin/plot_zvalue_featimps.R \\
        ${cluster_name}_phyloglm.csv ${feature} ${cluster_name}_zvalue_plot.png >> phyloglm.log 2>&1
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
    //fasta_ch = CLUSTERS_TO_FASTA().faa | flatten
    fasta_ch = Channel.fromPath("/vast/ri65fin/mcp_structures/nextflow/Caudoviricetes_Lake_Hot_unclustered/cluster_fastas/*.faa")

    // Call Mafft with fasta_ch, take fa output, flatten
    mafft_ch = MAFFT(fasta_ch).aligned | flatten

    // uncomment below if using a path that was previously generated
    // mafft_ch = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/Caudoviricetes_Lake_Hot_springs_split/mafft/*_mafft_aa.fa")

    //Call Fasttree with mafft_ch, take tree output, flatten
    fasttree_ch = FASTTREE(mafft_ch).tree | flatten
    // fasttree_ch = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/Caudoviricetes_Lake_Hot_springs_split/fasttree/*_combined_fastree_lg_pseudo.tree")

    // Optionally call TRIMAL, otherwise pass through mason_ch
    //msa_ch   = params.do_trimal ? TRIMAL(mafft_ch).faa : mafft_ch

    // Call RF_MODEL with mafft_ch
    msa_ch = params.do_trimal ? TRIMAL(mafft_ch).faa : mafft_ch
    rf_ch  = RF_MODEL(msa_ch)
    rf_ch_features = rf_ch.features | flatten
    rf_ch_data     = rf_ch.data_df  | flatten

    // rf_ch_features = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/Caudoviricetes_Lake_Hot_springs_split/RF_results/raw/RF_results/*_Feature_importance.tsv")
    // rf_ch_data = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/Caudoviricetes_Lake_Hot_springs_split/RF_results/raw/RF_results/*_data_df.tsv")

    //phylogml
    phyloglm_ch = PHYLOGLM(rf_ch_data, rf_ch_features, fasttree_ch)

    // Merge RF results
    MERGE_RF_RESULTS(
        rf_ch.features.collect(),
        rf_ch.reports.collect()
    )
}
