#!/usr/bin/env python

import argparse
import pandas as pd
import tqdm
import os
from Bio import SeqIO
from Bio.Phylo.PAML import yn00

#cl
#python script/dNdS.py -d input/replicons -T 23-01-25_results/TATracking_final_results.tsv --tblastn 23-01-25_results/1-tblastn/all_tblastn_raw_with_core.csv -o 23-01-25_results/dNdS


#script which aims to extract the coordinates of each TAs, to submit their to PAML yn00 soft to calculate the dN/dS
def main():
    
    args = get_args()
    df_TATracking = pd.read_csv(args.TATracking, sep = "\t", comment = "#")
    df_tblastn = pd.read_csv(args.tblastn, sep = "\t", comment = "#")
    prep_outdir(args.outdir)
    """"""
    if args.mode == "spot":
        d_features = prep_spot_extraction(df_TATracking, df_tblastn)
        df_TAs_features = pd.DataFrame.from_dict(d_features, orient = "index")
        df_TAs_features.to_csv(f"{args.outdir}/TAs_coordinates.tsv", sep = "\t")
        
        extract_seq(df_TAs_features, args.replicons, args.outdir)
    """"""
    
    yn00_handle = yn00.Yn00()
    yn00_handle.set_options(seqtype = 1)
        











def prep_spot_extraction(df_TATracking, df_tblastn):
    
    d_res = {}
    n = 0

    for TA in tqdm.tqdm(df_TATracking["Ref_TA"].drop_duplicates().tolist(), desc = "Extracting each TA within the same spot coordinates/features"):
        for genome in df_TATracking["Genome"].drop_duplicates().tolist():
            if df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["Same_core_spot"].tolist()[0] == "Yes":

                d_res[n] = {"TA": TA,
                            "genome": genome,
                            "contig": df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["TA_contig"].tolist()[0],
                            "Tox_leftc": df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["tox_left_coordinate"].tolist()[0],
                            "Tox_rightc": df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["tox_right_coordinate"].tolist()[0],
                            "Antitox_leftc": df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["antitox_left_coordinate"].tolist()[0],
                            "Antitox_rightc": df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["antitox_right_coordinate"].tolist()[0],
                            "strand": check_strand(TA, df_tblastn,
                                                   df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["TA_contig"].tolist()[0], 
                                                   df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["Core_family_left"].tolist()[0],
                                                   df_TATracking[(df_TATracking["Ref_TA"] == TA) & (df_TATracking["Genome"] == genome)]["Core_family_right"].tolist()[0])}
                
                n += 1


    return d_res






def check_strand(TA, df_tblastn, contig, left_core, right_core):
    
    tox_name = TA.split("-")[0]
    tmp_df = df_tblastn[(df_tblastn["qseqid"] == tox_name) & 
                        (df_tblastn["sseqid"] == contig) &
                        (df_tblastn["left_core_family"] == left_core) &
                        (df_tblastn["right_core_family"] == right_core)]
    
    if int(tmp_df["sstart"].tolist()[0]) < int(tmp_df["send"].tolist()[0]):
        return "+"
    else:
        return "-"
    


def extract_seq(df_TAs_features, replicon_folder, outdir):

    """
    Create two fna files for each TA, one with the toxins sequences, the other one with the antitoxin sequences
    """
    d_seqio = {}
    for genome in df_TAs_features["genome"].drop_duplicates().tolist():
        d_seqio[genome] = SeqIO.index(replicon_folder + "/"+ genome + ".fna", "fasta")
    
    for TA in tqdm.tqdm(df_TAs_features["TA"].drop_duplicates().tolist(), desc = "Extracting TAs sequences"):
        tox_str = ""
        antitox_str = ""
        for genome in df_TAs_features[df_TAs_features["TA"] == TA]["genome"].tolist():
            
            tox_seq = d_seqio[genome][df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["contig"].tolist()[0]].seq[
                int(df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["Tox_leftc"].tolist()[0])-1:
                int(df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["Tox_rightc"].tolist()[0])]

            antitox_seq = d_seqio[genome][df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["contig"].tolist()[0]].seq[
                int(df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["Antitox_leftc"].tolist()[0])-1:
                int(df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["Antitox_rightc"].tolist()[0])]
            
            if df_TAs_features[(df_TAs_features["TA"] == TA) & (df_TAs_features["genome"] == genome)]["strand"].tolist()[0] == "-":
                tox_seq = tox_seq.reverse_complement()
                antitox_seq = antitox_seq.reverse_complement()
            
            tox_str += f">{TA.split('-')[0]}_{genome}\n{tox_seq}\n"
            antitox_str += f">{TA.split('-')[1]}_{genome}\n{antitox_seq}\n"
        
        with open(f"{outdir}/TAs_sequences/{TA}_toxins_spot.fasta", "w") as f:
            f.write(tox_str)
            
        with open(f"{outdir}/TAs_sequences/{TA}_antitoxins_spot.fasta", "w") as f:
            f.write(antitox_str)
            






def prep_outdir(outdir):
    
    list_outfolders = [outdir, outdir + "/TAs_sequences", outdir + "/TAs_dNdS"]
    
    for i in list_outfolders:
        try :
            os.makedirs(i)
        except FileExistsError:
            continue
    
    


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--replicons",
                        help = "(REQUIRED) Directory all replicons files (fna format) used during the TATracking", required = True)
    parser.add_argument("-T", "--TATracking",
                        help = "TATracking result file", required = True)
    parser.add_argument("--tblastn",
                        help = "File generated by TATracking named all_tblastn_raw_core.csv", required = True)
    parser.add_argument("-o", "--outdir",
                        help = "Outdir to save the results", required = True)
    parser.add_argument("--mode",
                        help = "Depending on the kind of analysis you want: spot => calculate the dN/dS for each TAs only within its spot, all => cluster all TAs with mmseqs (default : spot)",
                        default = "spot")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()