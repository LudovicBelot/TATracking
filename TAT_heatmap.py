import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
from datetime import date
import argparse

#commandline : python script/TAT_heatmap.py -f 23-01-23_results/TATracking_final_results.tsv -o 23-01-23_results --order input/ordered_list_genomes_rooted_xbmm.lst -r PhAl.1022.00001 --names input/redo_name_dict.txt -t input/list_TA/TAlist_stabilizing.txt
# python script/TAT_heatmap.py -f 24-03-28_results_only_samecore/TATracking_final_results.tsv -o 24-03-28_results_only_samecore --order input/ordered_list_genomes_rooted_xbmm.lst -r PhAl.1022.00001 --names input/redo_name_dict.txt 

def main():

    plt.style.use("seaborn")

    args = get_args()
    df_TAT, outdir, reference_genome, list_genome, d_scientific_name, list_TA = format_args(args)

    d_res_TA_conserved_percentage = {} # index = strains genomes , columns = TA name
    d_res_noTA_conserved_percentage = {}
    d_res_TA_alone = {} #in case there are no others genes within the same corespot with the TA, we need a supplementar category

    d_res_TA_conserved_percentage, d_res_noTA_conserved_percentage, d_res_TA_alone = init_dict_df(list_genome, [d_res_TA_conserved_percentage, d_res_noTA_conserved_percentage, d_res_TA_alone])
    d_number_neighbours_ref = {} # k = TA name, v = numbers of genes used in this study and considered as neighbours of this TA
    #list_genome = df_TAT["Genome"].drop_duplicates().tolist()

    for TA in list_TA:
        #first we determine whether there are others genes (different of the TA) within the same core spot (in the ref genome), 
        # otherwise we cannot do the heatmap because there were no genes to track within the phylogeny

        if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == reference_genome)]["list_neighbours_TA"].values[0] == "[]":
            d_number_neighbours_ref[TA] = 0
            for genome in list_genome:
                if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Same_core_spot"].values[0] in ["Yes"] :#update to change "?" to no
                    d_res_TA_alone[genome][TA] = 1
                    d_res_noTA_conserved_percentage[genome][TA] = float('nan')
                    d_res_TA_conserved_percentage[genome][TA] = float('nan')

                else :
                    d_res_TA_alone[genome][TA] = float('nan')
                    d_res_noTA_conserved_percentage[genome][TA] = float('nan')
                    d_res_TA_conserved_percentage[genome][TA] =  float('nan')

        else :
            d_number_neighbours_ref[TA] = len([x for x in df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == reference_genome)]["list_neighbours_TA"].values[0].replace("[","").replace("]","").split(",") if x])
            for genome in list_genome:
                if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Same_core_spot"].values[0] in ["Yes"]: #update to change "?" to no
                    d_res_TA_conserved_percentage[genome][TA] = float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])
                    d_res_noTA_conserved_percentage[genome][TA] = float('nan')
                    d_res_TA_alone[genome][TA] = float('nan')

                elif df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Is_TA_homolog"].values[0] == "Yes" and df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Same_core_spot"].values[0] in ["No", "?"]:
                    if df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0] != "-":
                        d_res_noTA_conserved_percentage[genome][TA] = float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])
                        d_res_TA_conserved_percentage[genome][TA] = float('nan')
                        d_res_TA_alone[genome][TA] = float('nan')
                    else :
                        d_res_noTA_conserved_percentage[genome][TA] = float('nan')
                        d_res_TA_conserved_percentage[genome][TA] = float('nan')
                        d_res_TA_alone[genome][TA] = float('nan')

                elif df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["Is_TA_homolog"].values[0] in ["No"]: #update to change "?" to no
                    d_res_noTA_conserved_percentage[genome][TA] = float(df_TAT[(df_TAT["Ref_TA"] == TA) & (df_TAT["Genome"] == genome)]["%_conserved_genes/ref_genes"].values[0])
                    d_res_TA_conserved_percentage[genome][TA] = float('nan')
                    d_res_TA_alone[genome][TA] = float('nan')


    df_heatmap_neighbours_same_core_with_TA = pd.DataFrame.from_dict(d_res_TA_conserved_percentage, orient = "index").reindex(list_genome)
    df_heatmap_no_TA = pd.DataFrame.from_dict(d_res_noTA_conserved_percentage, orient = "index").reindex(list_genome)
    df_heatmap_TA_alone = pd.DataFrame.from_dict(d_res_TA_alone, orient = "index").reindex(list_genome)
    df_number_neighbours_ref = pd.DataFrame.from_dict(d_number_neighbours_ref, orient = "index")
    #df_number_neighbours_ref.to_csv(f"{outdir}/{date.today()}_n_neighbours_ref.tsv", sep ="\t")

    #calculate the number of TA conserved in each strain & the number of GI conserved ( based on the number of genes)
    d_data_per_genome = {}
    for genome in list_genome:
        d_data_per_genome[genome] = {"n_TA_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].count(), #+ df_heatmap_TA_alone.loc[genome].count(),
                                    "n_20%_neighbours_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].ge(20).sum()+df_heatmap_no_TA.loc[genome].ge(20).sum(),
                                    "n_40%_neighbours_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].ge(40).sum()+df_heatmap_no_TA.loc[genome].ge(40).sum(),
                                    "n_60%_neighbours_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].ge(60).sum()+df_heatmap_no_TA.loc[genome].ge(60).sum(),
                                    "n_80%_neighbours_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].ge(80).sum()+df_heatmap_no_TA.loc[genome].ge(80).sum(),
                                    "n_100%_neighbours_conserved": df_heatmap_neighbours_same_core_with_TA.loc[genome].ge(100).sum()+df_heatmap_no_TA.loc[genome].ge(100).sum()
                                    }

    df_data_per_genome = pd.DataFrame.from_dict(d_data_per_genome, orient = "index").reindex(list_genome)

    if d_scientific_name != None:
        df_heatmap_neighbours_same_core_with_TA.rename(index = d_scientific_name, inplace = True)
        df_heatmap_no_TA.rename(index = d_scientific_name, inplace = True)
        df_heatmap_TA_alone.rename(index = d_scientific_name, inplace = True)
        df_data_per_genome.rename(index = d_scientific_name, inplace = True)

    df_data_per_genome.to_csv(f"{outdir}/{date.today()}_TATracking_data_per_genome.tsv", sep = "\t")

    #creating the heatmap
    #Will maybe think of a way to autosize the fig
    fig, ax = plt.subplots(figsize = (24,12))
    cbar_ax1 = fig.add_axes([.915, .50, .018, .385])
    cbar_ax2 = fig.add_axes([.915, .10, .018, .385])

    fig.text(0.97,.682, "TA detected\n(Same core)", horizontalalignment='center', verticalalignment='center', fontsize = "x-large")
    fig.text(0.97,.282, "TA absent", horizontalalignment='center', verticalalignment='center', fontsize = "x-large")
    color_palette_with_TA = sns.color_palette("blend:#ffb5a1,#FF0000", as_cmap=True)
    color_palette_noTA = sns.color_palette("blend:#FFFFE0,#FFB81C", as_cmap=True)
    color_palette_TAalone = sns.color_palette("blend:#FFFFFF,#ff25fc", as_cmap=True)
    sns.heatmap(df_heatmap_neighbours_same_core_with_TA, cmap = color_palette_with_TA, vmin = 0, vmax = 100, linewidths = 0.20, ax = ax, cbar_ax= cbar_ax1 )
    sns.heatmap(df_heatmap_TA_alone, cmap = color_palette_TAalone, vmin = 0, vmax =1, linewidths = 0.20, cbar = False, ax = ax)
    sns.heatmap(df_heatmap_no_TA, cmap = color_palette_noTA, vmin = 0, vmax = 100, linewidths = 0.20, ax = ax, cbar_ax= cbar_ax2)
    plt.subplots_adjust(bottom=0.2, left = 0.175, top = 0.90)
    fig.savefig(f"{outdir}/{date.today()}_TATracking_heatmap.png")



def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", "-f",
                        help = " (REQUIRED) Input file generated from TATracking.py script", required = True)
    parser.add_argument("--reference", "-r",
                help = "(REQUIRED) Name of the reference genome used in this study (the one from which you provided the TA sequences)",
                required = True)
    parser.add_argument("--outdir", "-o",
                help ="(REQUIRED) Out directory in which the heatmap will be created",
                required = True)
    parser.add_argument("--order",
                    help = "(Optionnal) File with the genomes ordered (1 per line), this will be used to sort the rows of the heatmap in the descending order, default : None",
                    default = None)
    parser.add_argument("--names", "-n",
                    help = "(Optionnal) File with the 'informatic' names generated by PanACoTA and the scientific names separated by ':' (one per line), default : None",
                    default = None)
    parser.add_argument("--TAlist", "-t",
                        help = "File with the names of the TA to include in the heatmap (must be '{toxin_reference_name}-{antitoxin_reference_name}' (one per line), default = 'all'",
                        default = "all"
                        )

    args = parser.parse_args()
    return args


def format_args(args):
    # function which return a tuple with the needed information from the get_args function.

    df_TAT = pd.read_csv(args.file, index_col = 0, sep = "\t", comment = "#")
    outdir = args.outdir
    reference_genome = args.reference

    if args.order != None:
        with open(args.order, "r") as f:
            list_genome = [l.strip() for l in f if l]
    else :
        list_genome = df_TAT["Genome"].drop_duplicates().tolist()

    if args.names != None:
        d={}
        with open(args.names, "r") as f:
            for line in f:
                d[line.split(":")[0]] = line.split(":")[1].strip()
    else :
        d = None
    
    if args.TAlist != "all":
        TA_list = []
        with open(args.TAlist, "r") as f:
            for line in f:
                TA_list.append(line.strip())
    
    elif args.TAlist == "all":
        TA_list =  df_TAT["Ref_TA"].drop_duplicates().tolist()

    return df_TAT, outdir, reference_genome, list_genome, d, TA_list



def init_dict_df(list_genome, list_dict2init):
    
    res_d_list = []
    for d in list_dict2init:
        for genome in list_genome:
            d[genome] = {}
        res_d_list.append(d)

    return res_d_list




if __name__ == "__main__":
    main()