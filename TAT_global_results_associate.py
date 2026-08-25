import pandas as pd
import sys


#commandline 
# python script/TAT_global_results_associate.py 23-01-25_results/TATracking_final_results.tsv input/ordered_list_genomes_rooted_xbmm.lst input/redo_name_dict.txt 23-01-25_results/2023-06-20_results

#28-03-24 excluding '?'
#python script/TAT_global_results_associate.py 24-03-28_results_only_samecore/TATracking_final_results.csv input/ordered_list_genomes_rooted_xbmm.lst input/redo_name_dict.txt 24-03-28_results_only_samecore/4-global_results

def main():

    # small script which will create three files => a tsv file containing a dataframe of with for each TA (rows) whether or not the associated RGP is conserved or not
    # columns names are as follow {genome_name}_{%of RGP conservation} and one different column name "TA conserved"
    # cell values are either yes or no
    # second file is also a tsv file containing a dataframe which resume the data presented above as follows :
    # rows = genomes, columns = N_TAs with {X}% of their associated RGPs
    # finally last files will be :
    # rows = genomes, columns = N_TAs with {X}% of their associated RGPs / N_RGPs conserved with X% {with or without TAs}

    # NOTE: TAs immediatly surrounded by cores are excluded from this analyses and are reported in the excluded TA files 

    df_TATracking = pd.read_csv(sys.argv[1], sep = "\t")
    list_ordered_genomes = []
    with open(sys.argv[2], "r" ) as f:
        for line in f:
            list_ordered_genomes.append(line.strip())
    d_names = {}
    with open(sys.argv[3], "r" ) as f:
        for line in f:
            d_names[line.split(":")[0]] = line.split(":")[1].strip()

    outdir = sys.argv[4]



    list_TAs = df_TATracking[df_TATracking["%_conserved_genes/ref_genes"] != "-"]["Ref_TA"].drop_duplicates().tolist()
    list_TAs_excluded = df_TATracking[df_TATracking["%_conserved_genes/ref_genes"] == "-"]["Ref_TA"].drop_duplicates().tolist()
    #reporting which TAs are excluded 
    with open(f"{outdir}/TAs_surrounded by core.txt", "w") as f:
        str2write = "TAs excluded in this analysis:\n"
        for i in list_TAs_excluded:
            str2write += f"{i}\n"
        f.write(str2write)


    d_full_res = {}
    for TA in list_TAs:
        d_full_res[TA] = {}
        for genome in list_ordered_genomes:
            append_items2dict(check_RGP_conservation(TA, genome, df_TATracking), d_full_res[TA])

    # Counting also how many TAs surrounded by two core there are also !
    df_surroundedByCore_res = check_TAsurrounded(list_TAs_excluded, list_ordered_genomes, df_TATracking)
    df_surroundedByCore_res["genome_original_name"] = df_surroundedByCore_res.apply(lambda x: d_names[x["genome"]] , axis =1)
    df_surroundedByCore_res.to_csv(f"{outdir}/0_bis-TATracking_surroundedbycore_global_results.tsv", sep = "\t", index = False)

    df_full_res = pd.DataFrame.from_dict(d_full_res, orient = "index")
    df_full_res = df_full_res.applymap(lambda x: x[0])
    df_full_res.reset_index(inplace = True)
    df_full_res.rename(columns = {"index": "Ref_TA"}, inplace = True)
    df_full_res.to_csv(f"{outdir}/0-TATracking_full_res_associated_TA-RGP.tsv", sep = "\t", index = False)

    #working to produce the second file then:
    d2 = {}
    for genome in list_ordered_genomes:
        d2[genome] = {}
        for index_threshold in range(len([0,20,40,60,80,100])-1):
            d2[genome][f"[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%_ref_RGP"] = len(df_full_res[(df_full_res[f"TA_conserved_{genome}"] == "Yes") & (df_full_res[f"{genome}_RGP_[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%"] == "Yes")])


    df2_N_TA_RGP = pd.DataFrame.from_dict(d2, orient = "index")
    df2_N_TA_RGP.reset_index(inplace = True)
    df2_N_TA_RGP.rename(columns = {"index": "Genome"}, inplace = True)
    df2_N_TA_RGP["genome_original_name"] = df2_N_TA_RGP.apply(lambda x: d_names[x["Genome"]] , axis =1)
    df2_N_TA_RGP.to_csv(f"{outdir}/1-TATracking_N_TA_with_RGP.tsv", sep = "\t", index = False)

    #Now the last file with the ratio of TA_RGP_associated conserved divided by the number of RGP conserved total
    d3 = {}
    for genome in list_ordered_genomes:
        d3[genome] = {}
        for index_threshold in range(len([0,20,40,60,80,100])-1):
            try :
                d3[genome][f"Ratio_TA(RGP[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%)/n_RGP%"] = len(df_full_res[(df_full_res[f"TA_conserved_{genome}"] == "Yes") & (df_full_res[f"{genome}_RGP_[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%"] == "Yes")])/len(df_full_res[df_full_res[f"{genome}_RGP_[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%"] == "Yes"])
            except:
                d3[genome][f"Ratio_TA(RGP[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%)/n_RGP%"] = "Div0"

    df3= pd.DataFrame.from_dict(d3, orient = "index")
    df3.reset_index(inplace = True)
    df3.rename(columns = {"index": "Genome"}, inplace = True)
    df3["genome_original_name"] = df3.apply(lambda x: d_names[x["Genome"]] , axis =1)
    df3.to_csv(f"{outdir}/2-TATracking_ratio_TA-RGP_associated.tsv", sep = "\t", index = False)

    #Also creating a tsv file with the number of RGP conserved in each categorie for each TAs
    d_res_per_TA_all_RGPs = {}
    d_res_per_TA_onlywithTA_RGPs = {}

    for row in df_full_res.iterrows():
        d_res_per_TA_all_RGPs[row[0]] = {}
        d_res_per_TA_all_RGPs[row[0]]["Ref_TA"] = row[1]["Ref_TA"]
        d_res_per_TA_onlywithTA_RGPs[row[0]] = {}
        d_res_per_TA_onlywithTA_RGPs[row[0]]["Ref_TA"] = row[1]["Ref_TA"]

        count0_20_all = 0
        count21_40_all = 0
        count41_60_all = 0
        count61_80_all = 0
        count81_100_all = 0
        count0_20_withTA = 0
        count21_40_withTA = 0
        count41_60_withTA = 0
        count61_80_withTA = 0
        count81_100_withTA = 0

        for genome in list_ordered_genomes:
            if row[1][f"{genome}_RGP_[0-20]%"] == "Yes":
                count0_20_all += 1
                if row[1][f"TA_conserved_{genome}"] == "Yes":
                    count0_20_withTA += 1

            elif row[1][f"{genome}_RGP_[20-40]%"] == "Yes":
                count21_40_all += 1
                if row[1][f"TA_conserved_{genome}"] == "Yes":
                    count21_40_withTA += 1

            elif row[1][f"{genome}_RGP_[40-60]%"] == "Yes":
                count41_60_all += 1
                if row[1][f"TA_conserved_{genome}"] == "Yes":
                    count41_60_withTA += 1

            elif row[1][f"{genome}_RGP_[60-80]%"] == "Yes":
                count61_80_all += 1
                if row[1][f"TA_conserved_{genome}"] == "Yes":
                    count61_80_withTA += 1

            elif row[1][f"{genome}_RGP_[80-100]%"] == "Yes":
                count81_100_all += 1
                if row[1][f"TA_conserved_{genome}"] == "Yes":
                    count81_100_withTA += 1
        
        d_res_per_TA_all_RGPs[row[0]]["RGP_[0-20]%"] = count0_20_all
        d_res_per_TA_all_RGPs[row[0]]["RGP_[21-40]%"] = count21_40_all
        d_res_per_TA_all_RGPs[row[0]]["RGP_[41-60]%"] = count41_60_all
        d_res_per_TA_all_RGPs[row[0]]["RGP_[61-80]%"] = count61_80_all
        d_res_per_TA_all_RGPs[row[0]]["RGP_[81-100]%"] = count81_100_all

        d_res_per_TA_onlywithTA_RGPs[row[0]]["RGP_[0-20]%"] = count0_20_withTA
        d_res_per_TA_onlywithTA_RGPs[row[0]]["RGP_[21-40]%"] = count21_40_withTA
        d_res_per_TA_onlywithTA_RGPs[row[0]]["RGP_[41-60]%"] = count41_60_withTA
        d_res_per_TA_onlywithTA_RGPs[row[0]]["RGP_[61-80]%"] = count61_80_withTA
        d_res_per_TA_onlywithTA_RGPs[row[0]]["RGP_[81-100]%"] = count81_100_withTA


    df_res_per_TA_all_RGPs = pd.DataFrame.from_dict(d_res_per_TA_all_RGPs, orient = "index")
    df_res_per_TA_all_RGPs.to_csv(f"{outdir}/3-TATracking_global_results_per_TAs_all_RGPs.tsv", sep = "\t", index = False)

    df_res_per_TA_onlywithTA_RGPs = pd.DataFrame.from_dict(d_res_per_TA_onlywithTA_RGPs, orient = "index")
    df_res_per_TA_onlywithTA_RGPs.to_csv(f"{outdir}/4-TATracking_global_results_per_TAs_onlywithTAs_RGPs.tsv", sep = "\t", index = False)


def check_RGP_conservation(TA, genome, df_TATracking):

    d = {}
    df_tmp = df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]
    if df_tmp["Is_TA_homolog"].values[0] == "Yes":
        if df_tmp["Same_core_spot"].values[0] == "Yes" or df_tmp["Same_core_spot"].values[0] == "?":
            d[f"TA_conserved_{genome}"] = "Yes"
        else :
            d[f"TA_conserved_{genome}"] = "No"
    else :
        d[f"TA_conserved_{genome}"] = "No"

    for index_threshold in range (len([0,20,40,60,80,100])-1):
        if float(df_tmp["%_conserved_genes/ref_genes"].values[0]) > [0,20,40,60,80,100][index_threshold] and float(df_tmp["%_conserved_genes/ref_genes"].values[0]) <= [0,20,40,60,80,100][index_threshold+1]:
            d[f"{genome}_RGP_[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%"] = "Yes"        
        else :
            d[f"{genome}_RGP_[{[0,20,40,60,80,100][index_threshold]}-{[0,20,40,60,80,100][index_threshold+1]}]%"] = "No"
    
    if float(df_tmp["%_conserved_genes/ref_genes"].values[0]) == 0:
        d[f"{genome}_RGP_[0-20]%"] = "Yes"

    return d



def append_items2dict(d_in, d_out):
    for k, v in d_in.items():
        d_out.setdefault(k,[]).append(v)
    
    return d_out


def check_TAsurrounded(list_TAs, list_ordered_genomes, df_TATracking):
    d = {}
    for genome in list_ordered_genomes:
        d[genome] = {}
        for TA in list_TAs:
            df_tmp = df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]
            if df_tmp["Is_TA_homolog"].values[0] == "Yes":
                if df_tmp["Same_core_spot"].values[0] == "Yes" or df_tmp["Same_core_spot"].values[0] == "?":
                    d[genome][TA] = 1
                else :
                    d[genome][TA] = 0
            else :
                d[genome][TA] = 0
    
    df = pd.DataFrame.from_dict(d, orient = "index")
    list_columns = df.columns.tolist()
    df["total"] = df.sum(axis = 1)
    df["genome"] = df.index
    df = df[["genome" ,"total"]+list_columns]

    return df
    




if __name__ == "__main__":
    main()