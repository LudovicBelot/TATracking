import os
import sys
import pandas as pd
import matplotlib.pylab as plt

#commandline: python script/TAT_dist_TAcore.py input/core_genome/core_genome_photorhabdus_genus_features.csv input/gff3/PhAl.1022.00001.gff results/TATracking_final_results.tsv results/3-TAcore_distribution


#this script will return an histogramme representing the distribution of the number of genes within all corespot in the given genome and where the TA are located on.
#Note this script was created to work only on fully assembled genome

def main():

    #will maybe do a argparser later
    core_feature_file = sys.argv[1]
    gff_file = sys.argv[2]
    TA_tblastn_file = sys.argv[3]
    outdir = sys.argv[4]


    genome_name = gff_file.split("/")[-1].rsplit(".",1)[0]

    df_core_features = pd.read_csv(core_feature_file, sep = "\t", comment = "#", header = 0, index_col = 0)
    df_core_features = df_core_features[df_core_features["genome_name"] == genome_name].sort_values(by = "left_coordinate").reset_index(drop = False)
    df_tblastn_TA = pd.read_csv(TA_tblastn_file, sep = "\t", comment = "#", header = 0, index_col = 0)
    df_tblastn_TA = df_tblastn_TA[df_tblastn_TA["Genome"] == genome_name]
    for i in ["tox_left_coordinate","tox_right_coordinate","antitox_right_coordinate","antitox_left_coordinate"]:
        df_tblastn_TA[i] = pd.to_numeric(df_tblastn_TA[i])
    
    gff_columns = ["contig","source","type", "left_coordinate", "right_coordinate", "dot", "strand", "0" ,"description"]
    df_gff = pd.read_csv(gff_file, sep = "\t", comment = "#", header = 0, names = gff_columns)
    df_gff["description"] = df_gff["description"].apply(lambda x: x.split(";")[0].strip("ID="))
    TA_list = df_tblastn_TA["Ref_TA"].drop_duplicates().tolist()

    d_corespot = {} #k = number of the core spot v= {left_core: ..., right_core:..., numbers of genes within the core spot:... , number_TA:...}
    n_corespot = 1

    for core_gene_row in df_core_features.iterrows():
        if core_gene_row[0]+1 != len(df_core_features):
            left_core = core_gene_row[1]["index"]
            corespot_left_coordinate = core_gene_row[1]["left_coordinate"]
            right_core = df_core_features.iloc[core_gene_row[0]+1]["index"]
            corespot_right_coordinate = df_core_features.iloc[core_gene_row[0]+1]["right_coordinate"]
            list_genes_within_corespot = df_gff[(df_gff["type"] == "CDS") & (df_gff["left_coordinate"] > corespot_left_coordinate) & (df_gff["right_coordinate"] < corespot_right_coordinate)]["description"].tolist()
            list_TA_within_corespot = df_tblastn_TA[(df_tblastn_TA["tox_left_coordinate"] > corespot_left_coordinate) & (df_tblastn_TA["antitox_left_coordinate"] > corespot_left_coordinate) & (df_tblastn_TA["tox_right_coordinate"] < corespot_right_coordinate) & (df_tblastn_TA["tox_right_coordinate"] < corespot_right_coordinate)]["Ref_TA"].tolist()
            
        
        elif core_gene_row[0]+1 == len(df_core_features):
            left_core = core_gene_row[1]["index"]
            corespot_left_coordinate = core_gene_row[1]["right_coordinate"]
            right_core = df_core_features.iloc[0]["index"]
            corespot_right_coordinate = df_core_features.iloc[0]["left_coordinate"]
            list_genes_within_corespot = df_gff[(df_gff["type"] == "CDS") & ((df_gff["right_coordinate"] > corespot_left_coordinate) | (df_gff["left_coordinate"] < corespot_right_coordinate))]["description"].tolist()
            list_TA_within_corespot = df_tblastn_TA[((df_tblastn_TA["tox_right_coordinate"] < corespot_right_coordinate) & (df_tblastn_TA["antitox_right_coordinate"] < corespot_right_coordinate)) 
                                                    | ((df_tblastn_TA["tox_left_coordinate"] > corespot_left_coordinate) & (df_tblastn_TA["tox_left_coordinate"] > corespot_left_coordinate))]["Ref_TA"].tolist()



        d_corespot[n_corespot] = {"corespot_name": f"{left_core}-{right_core}", "left_coordinate": corespot_left_coordinate, "right_coordinate": corespot_right_coordinate, "n_TA": len(list_TA_within_corespot), "list_TA": (",").join(list_TA_within_corespot), "n_genes": len(list_genes_within_corespot), "list_genes": (",").join(list_genes_within_corespot)}
        n_corespot +=1  

    df_corespot = pd.DataFrame.from_dict(d_corespot, orient = "index")
    df_corespot.to_csv(f"{outdir}/{genome_name}_detailled_TAcorespot_distribution.tsv", sep = "\t")

    d_summary = {} # k =index, v = {number_genes_corespot:..., n_core_spots: ..., n_TA:...} 
    vmax_number_genes_corespot = df_corespot["n_genes"].max()
    
    for i in range(vmax_number_genes_corespot+1):
        d_summary[i] =  {"number_genes_within_corespot": i,
                        "n_core_spots": len(df_corespot[df_corespot["n_genes"] == i]),
                        "total_genes_in_these_corespot_lenght": len(df_corespot[df_corespot["n_genes"] == i])*i,
                        "n_TA": df_corespot[df_corespot["n_genes"] == i]["n_TA"].sum()
                        }

    df_summary_corespot = pd.DataFrame.from_dict(d_summary, orient = "index")
    df_summary_corespot.to_csv(f"{outdir}/{genome_name}_summary_TAcorespot_distribution.tsv", sep = "\t", index = False)




if __name__ == "__main__":
    main()