from plots import plot_report
from plots import plot_seq_logo
from plots import plot_AA_enrichment
from plots import plot_biome_feature_enrichment
from plots import plot_biome_features_full_protein
from plots import plot_protein_analysis

import pandas as pd
from pathlib import Path
import sys

if __name__ == '__main__':
    path_to_rf_results = Path(sys.argv[1])    #path to "RF_results" folder
    biome1 = sys.argv[2]
    biome2 = sys.argv[3]
    colour1 = sys.argv[4]
    colour2 = sys.argv[5]
    k = int(sys.argv[6])
    path_to_feat_imp = path_to_rf_results.joinpath("RandomForest_Report.tsv")
    df = pd.read_csv(path_to_feat_imp, sep='\t')

    # run the overall report plots
    plot_out_dir = path_to_rf_results.joinpath("Plots")
    Path.mkdir(plot_out_dir,exist_ok=True)

    # report_plot = plot_report.plot_report(df, biome1, biome2, colour1, colour2)
    # report_plot.savefig(plot_out_dir.joinpath("Rf_model_report_plots.png"))
    # report_plot.savefig(plot_out_dir.joinpath("Rf_model_report_plots.svg"))

    # p_v_r_plot = plot_report.plot_precision_vs_recall(df, biome1, biome2, colour1, colour2)
    # p_v_r_plot.savefig(plot_out_dir.joinpath("Rf_model_precision_v_recall.png"))

    for file in Path(path_to_rf_results.joinpath("raw")).glob("cluster_phyloglm_zvalue_new_01*_importance.tsv"):
        cluster_num = file.stem[0:-19]
        print(cluster_num)
        cluster_plot_dir = plot_out_dir.joinpath(f"{cluster_num}")
        Path.mkdir(cluster_plot_dir,exist_ok=True)

        plt = plot_seq_logo.plot(cluster_num, biome1, biome2, colour1, colour2,k)
        plt.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_seqLogo.png"))

        plt = plot_biome_feature_enrichment.plot(cluster_num, biome1, biome2, colour1, colour2,k)
        plt.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_biome_feat_enrich_bars.png"))

        plt = plot_biome_features_full_protein.plot(cluster_num, biome1, biome2, colour1, colour2,k)
        plt.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_biome_feat_enrich_bars_fullprot.png"))

        plt = plot_AA_enrichment.radar_biome_plot(cluster_num, colour1,colour2,top_k=k)
        plt.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_AA_radars.png"))

        plt = plot_AA_enrichment.AA_bar_plot(cluster_num, (biome1,biome2),top_k=k)
        plt.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_AA_enrich_bars.png"))

        # plt1 = plot_protein_analysis.main_plot(cluster_num, colour1, colour2,flex_window=9, polar_window=9,charge_window=20)
        # plt1, plt2, plt3, plt4 = plot_protein_analysis.main_plot(cluster_num, colour1, colour2,flex_window=9, polar_window=9,charge_window=20)
        # plt1.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_protein_analysis_subplots.png"), bbox_inches="tight")
        # plt1.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_protein_analysis_subplots.svg"), bbox_inches="tight")

        # plt2.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_protein_scale_flexibility.png"),bbox_inches="tight")
        # plt3.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_protein_scale_GRAVY.png"), bbox_inches="tight")
        # plt4.savefig(cluster_plot_dir.joinpath(f"{cluster_num}_protein_scale_chargeatpH7.png"), bbox_inches="tight")

        print(f"done... {cluster_num}")