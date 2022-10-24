import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

#commandline : python script/TAT_heatmap.py results/TATracking_final_results_only_0179-180.csv results input/ordered_list_genomes.lst


def main():
    #first, we need to create a dataframe with the results needed for the heatmap
    #will create a argparser later
    TAT_file = sys.argv[1]
    outdir = sys.argv[2]
    list_genome_file = sys.argv[3] #list file, one genome per line => will keep the order for the results
    plt.style.use("seaborn")

    #first need to format our results to get a dataframe with :
    #row : strain name
    #column TA name
    df_TAT = pd.read_csv(TAT_file, index_col = 0, sep = "\t", comment = "#")
    with open(list_genome_file, "r") as f:
        list_genome = [l.strip() for l in f if l]

    print(list_genome)
    d_res_TA_conserved_percentage = {} # index = strains genomes , columns = TA name
    d_res_noTA_conserved_percentage = {}
    #list_genome = df_TAT["Genome"].drop_duplicates().tolist()

    for TA in df_TAT["Ref_TA"].drop_duplicates().tolist():
        for genome in list_genome:
            if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Same_core_spot"].values[0] in ["Yes", "?"]:
                d_res_TA_conserved_percentage[genome] = {TA: float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])}
                d_res_noTA_conserved_percentage[genome] = {TA: None}

            elif df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Is_TA_homolog"].values[0] == "Yes" and df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Same_core_spot"].values[0] == "No":
                if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0] != "-":
                    d_res_noTA_conserved_percentage[genome] = {TA: float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])}
                    d_res_TA_conserved_percentage[genome] = {TA: None}
                else :
                    d_res_noTA_conserved_percentage[genome] = {TA: None}
                    d_res_TA_conserved_percentage[genome] = {TA: None}

            elif df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Is_TA_homolog"].values[0] == "No":
                d_res_noTA_conserved_percentage[genome] = {TA: float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])}
                d_res_TA_conserved_percentage[genome] = {TA: None}


    df_heatmap_neighbours_same_core_with_TA = pd.DataFrame.from_dict(d_res_TA_conserved_percentage, orient = "index").reindex(list_genome)
    print(df_heatmap_neighbours_same_core_with_TA)
    df_heatmap_no_TA = pd.DataFrame.from_dict(d_res_noTA_conserved_percentage, orient = "index").reindex(list_genome)
    print(df_heatmap_no_TA)

    #creating the heatmap
    fig, ax = plt.subplots(figsize = (12,12))
    color_palette_with_TA = sns.color_palette("blend:#ffb5a1,#FF0000", as_cmap=True)
    color_palette_noTA = sns.color_palette("blend:#FFFFE0,#FFB81C", as_cmap=True)
    sns.heatmap(df_heatmap_neighbours_same_core_with_TA, cmap = color_palette_with_TA, vmin = 0, vmax = 100, linewidths = 0.20)
    fig = sns.heatmap(df_heatmap_no_TA, cmap = color_palette_noTA, vmin = 0, vmax = 100, linewidths = 0.20).get_figure()
    fig.savefig(f"{outdir}/test.png")


if __name__ == "__main__":
    main()