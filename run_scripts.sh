# --- MODELS --- #
# for f in foldmason/*_db_out; do
#     clust_num=$(basename "$f" _db_out)
#     python3 bin/rf_model.py "${f}/${clust_num}_MSTA_aa.fa" /Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/amino_acid_analyis/imgvr_mgnify_combined/best_biome_per_mcp_protein.csv
# done

# --- PLOTTING --- #
# path_to_rf_results = Path(sys.argv[1])    #path to "RF_results" folder
# biome1 = sys.argv[2]
# biome2 = sys.argv[3]
# colour1 = sys.argv[4]
# colour2 = sys.argv[5]
# k = int(sys.argv[6])  # number of features

# python3 bin/plotting_script.py RF_results Lake Oceanic "#9dcc7e" "#2B8CBF" 100

# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300009526_000008|3300009526|Ga0115004_100234614" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300012282_000044|3300012282|Ga0157136_10002142" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "MGYP001198978755" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300020042_009324|3300020042|Ga0206640_10135673" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "MGYP002641143952" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "MGYP003704062755" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300043786_000143|3300043786|Ga0453106_0029613_1252_2892" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "MGYP003333885128" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300027697_000815|3300027697|Ga0209033_10097072" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300012758_000284|3300012758|Ga0138285_11887812" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm "IMGVR_UViG_3300029202_000060|3300029202|Ga0167843_1021547" RF_results

# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300009526_000008|3300009526|Ga0115004_100234614" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300012282_000044|3300012282|Ga0157136_10002142" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "MGYP001198978755" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300020042_009324|3300020042|Ga0206640_10135673" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "MGYP002641143952" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "MGYP003704062755" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300043786_000143|3300043786|Ga0453106_0029613_1252_2892" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "MGYP003333885128" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300027697_000815|3300027697|Ga0209033_10097072" RF_results
# python3 bin/map_seq_position.py cluster_0_1_mixed "IMGVR_UViG_3300012758_000284|3300012758|Ga0138285_11887812" RF_results

# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "EC6098_reference" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue_negative "EC6098_reference" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue_positive "EC6098_reference" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue_01 "EC6098_reference" RF_results

python3 bin/map_seq_position.py cluster_0_1_mixed_MSTA_aa_plus_ec6098re "MGYP003333944615" RF_results 
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue_new_01 "EC6098_reference" RF_results 


# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300009526_000008|3300009526|Ga0115004_100234614" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300012282_000044|3300012282|Ga0157136_10002142" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "MGYP001198978755" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300020042_009324|3300020042|Ga0206640_10135673" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "MGYP002641143952" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "MGYP003704062755" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300043786_000143|3300043786|Ga0453106_0029613_1252_2892" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "MGYP003333885128" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300027697_000815|3300027697|Ga0209033_10097072" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300012758_000284|3300012758|Ga0138285_11887812" RF_results
# python3 bin/map_seq_position.py cluster_phyloglm_zvalue "IMGVR_UViG_3300029202_000060|3300029202|Ga0167843_1021547" RF_results