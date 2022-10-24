from collections import OrderedDict
from weakref import ref
import pandas as pd
from datetime import date
import sys
import os
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
from Bio.Blast.Applications import NcbitblastnCommandline, NcbiblastpCommandline


#commandline test : python script/TATracking.py input/TAseq.fna input/replicons input/core_genome/core_genome_photorhabdus_genus_features.csv PhAl.1022.00001 input


def main():

    #will probably change the sys method to argparse later
    TA_fasta_file = sys.argv[1]
    genome_folder = os.path.expanduser(sys.argv[2])
    core_features_file = sys.argv[3]
    reference_genome = sys.argv[4]
    protein_gff_folder = os.path.expanduser(sys.argv[5])

    od_TA = TAfna2dic(TA_fasta_file)
    TAfile_seqIO = SeqIO.index(TA_fasta_file, "fasta")
    list_genomes = []

    #create_outdir()
    tmp_folder = "tmp"
    outdir = "results"

    #First create (if not already done) a nucleotidic blast db for each genome to study
    for file in os.listdir(genome_folder):
        list_genomes.append(file) #useful later to perform the tblastn analysis in each of these genomes
        if (file.endswith(".fna") or file.endswith(".fasta")) and os.path.exists(tmp_folder+'/blastdb/'+file+".nhr") == False:
            os.system(f"makeblastdb -in {genome_folder+'/'+file} -out {tmp_folder}/blastdb/{file} -parse_seqids -dbtype nucl")

    #now for each of these genomes we will search for the TA using tblastn
    print("Performing tblastN for each TAs")
    for TA_operon in tqdm(od_TA.values()):
        TA_tblastn(TA_operon, list_genomes, TAfile_seqIO, tmp_folder, f"{outdir}/1-tblastn")

    #create a tsv file containing all tblastn from the analysis + their localization compared to the core genome (for each hit)
    print("Associating each toxin/antitoxin hits with their corespot location")
    localize_with_core(outdir, core_features_file) #uncomment here

    #Then associate each Toxin hit with an antitoxin hit if within a same core spot and with an intergenic interval max of 150 bp (?)
    #Note from now, we keep only hits (toxin & antitoxin) with a %id >= 80 and a %cov >= 80%
    print("Determining the localization of each complete TA hits")
    df_full_TA_tblastn = associate_TA_tblastn_hits(f"{outdir}/1-tblastn/all_tblastn_raw_with_core.csv", od_TA) #uncomment here
    df_full_TA_tblastn.to_csv(f"{outdir}/1-tblastn/full_TA_tblastn_with_core.csv" , sep = "\t")

    # for each TA operon hit, we get the neigbouring_genes genes in each genome
    #df_ref_neigbouring_genes = get_ref_neighbouring_genes(df_full_TA_tblastn, reference_genome, protein_gff_folder)
    #df_ref_neigbouring_genes.to_csv(f"{outdir}/2-neighbouring_genes/ref_genome_neigbouring_genes.csv", sep = "\t")
    print("Getting the neighbouring genes of each TA hits")
    df_all_neigbouring_genes = get_all_neigbouring_genes(df_full_TA_tblastn, protein_gff_folder)
    df_all_neigbouring_genes.to_csv(f"{outdir}/2-neighbouring_genes/all_genome_neigbouring_genes.csv", sep = "\t")

    # now we do blast the TAs neighbouring genes from the reference genome in with all neighbouring genes 
    print("Blastp of the TAs neighbouring genes with the reference neighbouring genes")
    create_db_blastp(df_all_neigbouring_genes, reference_genome, list_genomes, protein_gff_folder, tmp_folder)

    #for each TA, we do a blastP of each neigbouring genes from the reference genome with the neighbouring genes from the other genomes
    blastp_neigbouring_genes(df_all_neigbouring_genes, reference_genome, tmp_folder, outdir)
    #Then we also look for the neigbouring genes in the genome for which we do not have the TA within the same core spot
    print("Blastp of the genes within the same core spot (but without the TA) in each genome compared to the reference genome")
    d_corespot_testedgenes_noTA, d_genome_with_TA_same_core_spot = blastp_neigbouring_genes_noTA_genomes(df_all_neigbouring_genes, core_features_file, reference_genome, list_genomes, protein_gff_folder, tmp_folder, outdir)

    #analyzing the results of the previous steps
    print("Computing the results file")
    df_res = analyze_blastp_neighbouring_genes_step(df_all_neigbouring_genes, d_corespot_testedgenes_noTA, reference_genome, list_genomes, outdir)
    df_res.to_csv(f"{outdir}/TATracking_final_results.csv",sep ="\t")

############################################################################################################################################################


def TAfna2dic(fasta_file):
    od = OrderedDict()
    
    with open(fasta_file, "r") as f :
        count = 0
        n_TA = 1
        od[n_TA] = []
        for line in f:
            if line.startswith(">"):
                count +=1
                od[n_TA].append(line.replace(">","").replace("\n",""))
            if count == 2 :
                count =0
                n_TA +=1 
                od[n_TA] = []
    
    to_delete = []
    for k,v in od.items():
        if v == []:
            to_delete.append(k)
    for i in to_delete:
        od.pop(k)

    return od

def TA_tblastn(TA_operon, list_genomes, TAfile_seqIO, tmp_folder, outdir):
    # first need to create a tmp file with the seq of the toxin
    with open(tmp_folder+"/"+"tmp_tblastn_query.fasta","w") as f:
        f.write(f">{TAfile_seqIO[TA_operon[0]].id}\n{TAfile_seqIO[TA_operon[0]].seq}")

    str_columns_tblastn = "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    #now doing the tblastn 
    for genome in list_genomes :
        tblastn_cline = NcbitblastnCommandline(query = tmp_folder+"/"+"tmp_tblastn_query.fasta", db =tmp_folder+"/blastdb/"+genome, evalue = 0.1,
                                            outfmt = str_columns_tblastn , out= f"{outdir}/raw_data/{TA_operon[0]}_{genome.rsplit('.',1)[0]}.csv")
        stdout, stderr = tblastn_cline()

    # doing the same for the antitoxin 
    with open(tmp_folder+"/"+"tmp_tblastn_query.fasta","w") as f:
        f.write(f">{TAfile_seqIO[TA_operon[1]].id}\n{TAfile_seqIO[TA_operon[1]].seq}")

    #now doing the tblastn 
    for genome in list_genomes :
        tblastn_cline = NcbitblastnCommandline(query = tmp_folder+"/"+"tmp_tblastn_query.fasta", db =tmp_folder+"/blastdb/"+genome, evalue = 0.1,
                                            outfmt = str_columns_tblastn , out= f"{outdir}/raw_data/{TA_operon[1]}_{genome.rsplit('.',1)[0]}.csv")
        stdout, stderr = tblastn_cline()

def localize_with_core(outdir, core_file):

    #first getting the features of each core genes into a dataframe (from the provided core file)
    df_core = pd.read_csv(core_file, sep = "\t", comment = "#")

    #now looking for each hit of the TA tblastn and their localization compared to the core genome
    list_tblastn_df = []
    list_columns = ["qseqid", "qlen", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    for file in os.listdir(f"{outdir}/1-tblastn/raw_data"):
        df = pd.read_csv(f"{outdir}/1-tblastn/raw_data/{file}", sep = "\s+", names = list_columns)
        df_updated = df.copy(deep = True)

        #first we limit the core_dataframe to the genome & contig on which we found each tblastn hit
        for tblastn_hit in df.iterrows():
            df_updated.loc[tblastn_hit[0],"left_core"], df_updated.loc[tblastn_hit[0], "left_core_family"], df_updated.loc[tblastn_hit[0], "right_core"], df_updated.loc[tblastn_hit[0], "right_core_family"] = get_closest_core(tblastn_hit, df_core)
        list_tblastn_df.append(df_updated)    
    
    all_df_tblastn_raw_with_core = pd.concat(list_tblastn_df, ignore_index= True)
    all_df_tblastn_raw_with_core.to_csv(f"{outdir}/1-tblastn/all_tblastn_raw_with_core.csv", sep = "\t")

    return None

def get_closest_core(row_tblastn_hit, df_core):

    left_core_df = df_core[(df_core["contig"] == row_tblastn_hit[1]["sseqid"]) & (df_core["left_coordinate"] <= min(row_tblastn_hit[1]["sstart"], row_tblastn_hit[1]["send"]))].sort_values("left_coordinate")
    if left_core_df.empty == False:
        left_core_gene = left_core_df.iloc[-1]["Unnamed: 0"]
        left_core_family = int(left_core_df.iloc[-1]["core_family"])
    else :
        left_core_gene = "None"
        left_core_family = "None"

    right_core_df = df_core[(df_core["contig"] == row_tblastn_hit[1]["sseqid"]) & (df_core["right_coordinate"] >= max(row_tblastn_hit[1]["sstart"], row_tblastn_hit[1]["send"]))].sort_values("right_coordinate")
    if right_core_df.empty == False:
        right_core_gene = right_core_df.iloc[0]["Unnamed: 0"]
        right_core_family = int(right_core_df.iloc[0]["core_family"])
    else :
        right_core_gene = "None"
        right_core_family = "None"

    return left_core_gene, left_core_family, right_core_gene, right_core_family

def associate_TA_tblastn_hits(all_tblastn_csv , od_TA):
    
    #return a new dataframe with toxin/antitoxin tblastn hits associated if within the same core spot & not farthest than 150 bp
    df_all_tblastn = pd.read_csv(all_tblastn_csv, sep = "\t", index_col = 0)
    df_all_tblastn_80 = df_all_tblastn[(df_all_tblastn["pident"] >= 80) & (df_all_tblastn["length"]/df_all_tblastn["qlen"]*100 >= 80)]
    d_full_TA = {}
    n_index = 0

    for TA_operon in od_TA.values():
        for toxin_hit_row in df_all_tblastn_80[df_all_tblastn_80["qseqid"] == TA_operon[0]].iterrows():
            df_tmp_antitox = df_all_tblastn_80[(df_all_tblastn_80["qseqid"] == TA_operon[1]) & (df_all_tblastn_80["left_core_family"] == toxin_hit_row[1]["left_core_family"])
            & (df_all_tblastn_80["right_core_family"] == toxin_hit_row[1]["right_core_family"])]
            #now we keep only close hits >150bp
            if df_tmp_antitox.empty == False :
                df_tmp_antitox = df_tmp_antitox[(abs(df_tmp_antitox["sstart"] - min(int(toxin_hit_row[1]["sstart"]),int(toxin_hit_row[1]["send"]))) <=150) |
                                                (abs(df_tmp_antitox["sstart"] - max(int(toxin_hit_row[1]["sstart"]),int(toxin_hit_row[1]["send"]))) <=150) |
                                                (abs(df_tmp_antitox["send"] - min(int(toxin_hit_row[1]["sstart"]),int(toxin_hit_row[1]["send"]))) <=150) |
                                                (abs(df_tmp_antitox["send"] - max(int(toxin_hit_row[1]["sstart"]),int(toxin_hit_row[1]["send"]))) <=150) ]
                if df_tmp_antitox.empty == False:
                    #if there are still multiples hits, we keep the one we the highest e-value (I don't think it should happen often)
                    antitox_hit_row = df_tmp_antitox.sort_values("evalue").iloc[0]
                    d_full_TA[n_index] = {"ref_toxin": toxin_hit_row[1]["qseqid"], "tox_pident": toxin_hit_row[1]["pident"], "tox_pcov": toxin_hit_row[1]["length"]/toxin_hit_row[1]["qlen"]*100, "tox_evalue": toxin_hit_row[1]["evalue"],
                                        "tox_left_coordinate": min(toxin_hit_row[1]["sstart"], toxin_hit_row[1]["send"]), "tox_right_coordinate": max(toxin_hit_row[1]["sstart"],toxin_hit_row[1]["send"]),
                                        "ref_antitoxin": antitox_hit_row["qseqid"], "antitox_pident": antitox_hit_row["pident"], "antitox_pcov": antitox_hit_row["length"]/antitox_hit_row["qlen"]*100, "antitox_evalue": antitox_hit_row["evalue"],
                                        "antitox_left_coordinate": min(antitox_hit_row["sstart"], antitox_hit_row["send"]), "antitox_right_coordinate": max(antitox_hit_row["sstart"], antitox_hit_row["send"]),
                                        "genome_name": toxin_hit_row[1]["sseqid"][:-5], "sseqid_TA": toxin_hit_row[1]["sseqid"], "left_core_gene": toxin_hit_row[1]["left_core"], "left_core_family": toxin_hit_row[1]["left_core_family"],
                                        "right_core_gene": toxin_hit_row[1]["right_core"], "right_core_family": toxin_hit_row[1]["right_core_family"]
                                        }
                    n_index += 1
    
    df_full_TA = pd.DataFrame.from_dict(d_full_TA, orient = "index")
    return df_full_TA

def get_ref_neighbouring_genes(df_full_TA_tblastn, reference_genome, protein_gff_folder):
    #return a dataframe with the 10 closest genes (5 upstream & downstream) of the TA_tblastn coordinates only in the reference genome
    #Note: it will stop if one of these genes is part of the core genome 
    #Also: overlapping genes within the coordinates of the TA tblastn hit are excluded.

    gff_columns = ["contig","source","type", "left_coordinate", "right_coordinate", "dot", "strand", "0" ,"description"]

    df_full_TA_ref = df_full_TA_tblastn[df_full_TA_tblastn["genome_name"] == reference_genome]
    df_ref_genome = pd.read_csv(f"{protein_gff_folder}/gff3/{reference_genome}.gff", sep ="\t", names = gff_columns, comment ="#")
    df_ref_genome["description"] = df_ref_genome["description"].apply(lambda x: x.split(';')[0].split('=')[-1])
    
    d_TA_ref_neighbouring_genes = {}
    for TA_detected in df_full_TA_ref.iterrows():
        tmp_index_core_left = df_ref_genome[df_ref_genome["description"] == TA_detected[1]["left_core_gene"]].index[0]
        tmp_index_core_right = df_ref_genome[df_ref_genome["description"] == TA_detected[1]["right_core_gene"]].index[0]
        df_core_spot = df_ref_genome.iloc[tmp_index_core_left +1 :tmp_index_core_right]
        df_up_core_spot = df_core_spot[df_core_spot["right_coordinate"] <= min(TA_detected[1]["tox_left_coordinate"],TA_detected[1]["antitox_left_coordinate"])]
        df_dw_core_spot = df_core_spot[df_core_spot["left_coordinate"] >= max(TA_detected[1]["tox_right_coordinate"],TA_detected[1]["antitox_right_coordinate"])]

        list_ref_neighbouring_genes = []
        if len(df_up_core_spot) >= 5 :
            list_ref_neighbouring_genes = df_up_core_spot.iloc[-5:]["description"].tolist()
        elif len(df_up_core_spot) < 5 and len(df_up_core_spot) > 0: #In case there are not 5 genes upstream within the same core spot, we take the maximum we can 
            list_ref_neighbouring_genes = df_up_core_spot.iloc[-len(df_up_core_spot):]["description"].tolist()
        
        if len(df_dw_core_spot) >= 5 :
            list_ref_neighbouring_genes += df_dw_core_spot.iloc[:5]["description"].tolist()
        if len(df_dw_core_spot) < 5 and len(df_dw_core_spot) > 0:
            list_ref_neighbouring_genes += df_dw_core_spot.iloc[:len(df_dw_core_spot)]["description"].tolist()
        
        d_TA_ref_neighbouring_genes[f"{TA_detected[1]['ref_toxin']}-{TA_detected[1]['ref_antitoxin']}"] = {"neighbouring_genes" : ",".join(list_ref_neighbouring_genes),
                                                                                                            "n_neigbours": len(list_ref_neighbouring_genes)
                                                                                                            }
    
    df_TA_ref_neighbouring_genes = pd.DataFrame.from_dict(d_TA_ref_neighbouring_genes, orient = "index")
    
    return df_TA_ref_neighbouring_genes


def get_all_neigbouring_genes(df_full_TA_tblastn, protein_gff_folder):
    #return a dataframe with the 10 closest genes (5 upstream & downstream) of the TA_tblastn coordinates only in the reference genome
    #Note: it will stop if one of these genes is part of the core genome 
    #Also: overlapping genes within the coordinates of the TA tblastn hit are excluded.

    gff_columns = ["contig","source","type", "left_coordinate", "right_coordinate", "dot", "strand", "0" ,"description"]
    list_genomes = df_full_TA_tblastn["genome_name"].drop_duplicates().tolist()
    n_index = 0
    d_TA_hits_neighbours = {}

    for genome in list_genomes :
        df_gff_this_genome = pd.read_csv(f"{protein_gff_folder}/gff3/{genome}.gff", sep = "\t", comment = "#", names = gff_columns)
        df_gff_this_genome["description"] = df_gff_this_genome["description"].apply(lambda x: x.split(';')[0].split('=')[-1])
        df_gff_this_genome = df_gff_this_genome[df_gff_this_genome["type"] == "CDS"]
        df_TA_this_genome = df_full_TA_tblastn[df_full_TA_tblastn["genome_name"] == genome]

        for TA_hit in df_TA_this_genome.iterrows():
            if TA_hit[1]["left_core_gene"] != "None":
                tmp_index_core_left = df_gff_this_genome[df_gff_this_genome["description"] == TA_hit[1]["left_core_gene"]].index[0]
            else :
                tmp_index_core_left = None
            if TA_hit[1]["right_core_gene"] != "None":
                tmp_index_core_right = df_gff_this_genome[df_gff_this_genome["description"] == TA_hit[1]["right_core_gene"]].index[0]
            else :
                tmp_index_core_right = None
            
            if tmp_index_core_right != None and tmp_index_core_left != None :
                df_core_spot = df_gff_this_genome[df_gff_this_genome["contig"] == TA_hit[1]["sseqid_TA"]].loc[tmp_index_core_left +1 :tmp_index_core_right-1]
            elif tmp_index_core_right == None and tmp_index_core_left != None :
                df_core_spot = df_gff_this_genome[df_gff_this_genome["contig"] == TA_hit[1]["sseqid_TA"]].loc[tmp_index_core_left +1 :]
            elif tmp_index_core_right != None and tmp_index_core_left == None :
                df_core_spot = df_gff_this_genome[df_gff_this_genome["contig"] == TA_hit[1]["sseqid_TA"]].loc[:tmp_index_core_right-1]
            else :
                df_core_spot = df_gff_this_genome[df_gff_this_genome["contig"] == TA_hit[1]["sseqid_TA"]]


            list_neighbouring_genes = []
            if df_core_spot.empty == False: #in case it remains only the tblastn hits within  a core spot but no complete ORFs
                df_up_core_spot = df_core_spot[df_core_spot["right_coordinate"] <= min(TA_hit[1]["tox_left_coordinate"],TA_hit[1]["antitox_left_coordinate"])]
                df_dw_core_spot = df_core_spot[df_core_spot["left_coordinate"] >= max(TA_hit[1]["tox_right_coordinate"],TA_hit[1]["antitox_right_coordinate"])]

                if len(df_up_core_spot) >= 5 :
                    list_neighbouring_genes = df_up_core_spot.iloc[-5:]["description"].tolist()
                elif len(df_up_core_spot) < 5 and len(df_up_core_spot) > 0: #In case there are not 5 genes upstream within the same core spot, we take the maximum we can 
                    list_neighbouring_genes = df_up_core_spot.iloc[-len(df_up_core_spot):]["description"].tolist()
                
                if len(df_dw_core_spot) >= 5 :
                    list_neighbouring_genes += df_dw_core_spot.iloc[:5]["description"].tolist()
                if len(df_dw_core_spot) < 5 and len(df_dw_core_spot) > 0:
                    list_neighbouring_genes += df_dw_core_spot.iloc[:len(df_dw_core_spot)]["description"].tolist()

            d_TA_hits_neighbours[n_index] = {"TA_homolog_of": f"{TA_hit[1]['ref_toxin']}-{TA_hit[1]['ref_antitoxin']}",
                                            "genome_name": genome,
                                            "contig": TA_hit[1]["sseqid_TA"],
                                            "left_core_gene": TA_hit[1]["left_core_gene"],
                                            "left_core_family": TA_hit[1]["left_core_family"],
                                            "right_core_gene": TA_hit[1]["right_core_gene"],
                                            "right_core_family": TA_hit[1]["right_core_family"],
                                            "n_neighbours": len(list_neighbouring_genes),
                                            'neighbours_genes': ",".join(list_neighbouring_genes)
                                            }
            n_index +=1

    df_TA_hits_neighbours = pd.DataFrame.from_dict(d_TA_hits_neighbours, orient = "index")
    #df_TA_hits_neighbours.sort_values(by = "TA_homolog_of", inplace = True)
    return df_TA_hits_neighbours

def create_db_blastp(df_all_neigbouring_genes, reference_genome, list_genomes, protein_gff_folder, tmp_folder): 

    #create two .faa files, one with all neighbouring genes from the the reference genome, and one with all neighbouring genes within the others genomes
    df_ref = df_all_neigbouring_genes[df_all_neigbouring_genes["genome_name"] == reference_genome]
    str_list_ref_neighbours = ",".join(df_ref["neighbours_genes"].tolist())
    list_ref_neighbours =  [x for x in str_list_ref_neighbours.split(",") if x] #create a list without empty string
    list_ref_neighbours = list(set(list_ref_neighbours)) #delete duplicates

    ref_seqIO = SeqIO.index(f"{protein_gff_folder}/Proteins/{reference_genome}.prt", "fasta")
    with open(f"{tmp_folder}/ref_genome_neighbours_db.faa", "w") as f:
        for neighbor in list_ref_neighbours:
            f.write(f">{ref_seqIO[neighbor].id}\n{ref_seqIO[neighbor].seq}\n")

    #os.system(f"makeblastdb -in {tmp_folder}/ref_genome_db.faa -out {tmp_folder}/ref_genome_db -parse_seqids -dbtype prot")

    """
    #now creating the database of the neighbouring genes in the others genomes
    df_others = df_all_neigbouring_genes[df_all_neigbouring_genes["genome_name"] != reference_genome]
    str_list_others_neighbours =",".join(df_others["neighbours_genes"].tolist())
    list_others_neighbours = [x for x in str_list_others_neighbours.split(",") if x]
    list_others_neighbours =list(set(list_others_neighbours))
    list_others_genomes = df_all_neigbouring_genes[df_all_neigbouring_genes["genome_name"] != reference_genome]["genome_name"].drop_duplicates().tolist()

    all_seqIO = None
    for g in list_others_genomes:
        if all_seqIO == None :
            all_seqIO = SeqIO.index(f"{protein_gff_folder}/Proteins/{g}.prt", "fasta")
        else :
            all_seqIO = MultiIndexDict(all_seqIO, SeqIO.index(f"{protein_gff_folder}/Proteins/{g}.prt", "fasta"))
    
    with open(f"{tmp_folder}/all_others_genome_neighbours_db.faa", "w") as f:
        for neighbor in list_others_neighbours:
            f.write(f">{all_seqIO[neighbor].id}\n{all_seqIO[neighbor].seq}\n")

    os.system(f"makeblastdb -in {tmp_folder}/all_others_genome_neighbours_db.faa -out {tmp_folder}/all_others_genome_neighbours_db -parse_seqids -dbtype prot")
    """
    #now creating a prt database with all proteomes sequences except for the reference genome
    list_genomes = [x.rsplit(".",1)[0] for x in list_genomes if x.rsplit(".",1)[0] != reference_genome]
    str_cat_proteome = ""
    for g in list_genomes:
        str_cat_proteome += f"{protein_gff_folder}/Proteins/{g}.prt "
    
    os.system(f"cat {str_cat_proteome} > {tmp_folder}/all_proteome_db.prt")
    os.system(f"makeblastdb -in {tmp_folder}/all_proteome_db.prt -out {tmp_folder}/all_proteome_db -parse_seqids -dbtype prot")



def blastp_neigbouring_genes(df_all_neigbouring_genes, reference_genome, tmp_folder, outdir):

    #recreating the SeqIO index iterator with all genes neighbouring the TAs in the 
    ref_neigbours_seqio = SeqIO.index(f"{tmp_folder}/ref_genome_neighbours_db.faa", "fasta")

    #doing the blastp analysis to compare each neigbouring genes from the reference genome to the ones of the other genomes
    list_columns = ["qseqid", "qlen", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    str_columns_blastp = "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    for TA in df_all_neigbouring_genes["TA_homolog_of"].drop_duplicates().tolist():
        #first generate a temporary file with the sequence of the TA neigbouring proteins
        list_ref_neigbours_this_TA = df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == reference_genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)].drop_duplicates(subset = "TA_homolog_of")["neighbours_genes"].values[0]
        list_ref_neigbours_this_TA = [x for x in list_ref_neigbours_this_TA.split(",") if x]
        with open(f"{tmp_folder}/TA_neigbouring_genes_tmp_blastp_query.faa", "w") as tmp_query:
            for i in list_ref_neigbours_this_TA:
                tmp_query.write(f">{ref_neigbours_seqio[i].id}\n{ref_neigbours_seqio[i].seq}\n")
        
        #now doing the blastp with these protein sequences and restricting the blast subject database to all proteins sequences which are neigbouring the TA in the others genome.
        #for restricting the search it needs a file with the id of all sequences to use, one per line
        list_other_neigbours_this_TA = df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] != reference_genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]["neighbours_genes"].tolist()
        list_other_neigbours_this_TA_formated = []
        for list_neighbours in list_other_neigbours_this_TA:
            list_other_neigbours_this_TA_formated += [x for x in list_neighbours.split(",") if x]
        
        with open(f"{tmp_folder}/tmp_seqidlist.lst", "w") as tmp_seqidlist:
            for i in list_other_neigbours_this_TA_formated:
                tmp_seqidlist.write(f"{i}\n")


        blastp_cline = NcbiblastpCommandline(query = f"{tmp_folder}/TA_neigbouring_genes_tmp_blastp_query.faa", db = f"{tmp_folder}/all_proteome_db", evalue = 0.1,
                                            outfmt = str_columns_blastp , seqidlist = f"{tmp_folder}/tmp_seqidlist.lst",
                                            out= f"{outdir}/2-neighbouring_genes/blastp_raw_data/{TA}_blastp_neighbours_near_TA.csv")
        stdout, stderr = blastp_cline()


def blastp_neigbouring_genes_noTA_genomes(df_all_neigbouring_genes, core_features_file, reference_genome, list_genomes, protein_gff_folder, tmp_folder, outdir):

    list_genomes = [x.rsplit(".",1)[0] for x in list_genomes if x.rsplit(".",1)[0] != reference_genome]

    #First we store the core spot of each TA in the ref genome in a dictionnary k = TA => v: [core_gene_family_left, core_gene_family_right]
    d_TA_core_ref = {}
    for TA in df_all_neigbouring_genes.drop_duplicates(subset = "TA_homolog_of")["TA_homolog_of"].tolist():
        d_TA_core_ref[TA] = [df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == reference_genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)].drop_duplicates(subset = "TA_homolog_of")["left_core_family"].values[0],
                            df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == reference_genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)].drop_duplicates(subset = "TA_homolog_of")["right_core_family"].values[0]]

    #now checking if the TA found in the others genomes are part of the same core spot (if so, we don't need to do the blastp analysis here for those genomes)
    d_genome_with_TA_same_core_spot = {} # k = TA, v = list of genomes with the TA within the same core spot 
    for TA in d_TA_core_ref.keys():
        d_genome_with_TA_same_core_spot[TA] = []
        for genome in list_genomes :
            df_tmp = df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]
            for row in df_tmp.iterrows():
                if (row[1]["left_core_family"] in d_TA_core_ref[TA] and row[1]["right_core_family"] in d_TA_core_ref[TA]) or (row[1]["left_core_family"] == "None" and row[1]["right_core_family"] in d_TA_core_ref[TA]) or (row[1]["left_core_family"] in d_TA_core_ref[TA] and row[1]["right_core_family"] == "None"):
                    d_genome_with_TA_same_core_spot[TA].append(genome)
                    break
    
    #creating a dictionnary with the genomes for which we do not have the TA within the same core spot
    d_genome_with_noTA_same_core_spot = {}
    for k,v in d_genome_with_TA_same_core_spot.items():
        d_genome_with_noTA_same_core_spot[k] = [x for x in list_genomes if x not in v]

    # And finally getting the genes within the same core spot in the genomes stored dict "d_genome_with_noTA_same_core_spot"
    # Note: in case the core spot is incomplete, meaning that both core spot are located on different contigs, 
    # we're trying to determine first if they are other core genes on the same contig which could help us to determine the core spot,
    # if it's not possible we're trying the 4 combinaisons and keeping the best results

    #storing in a dict all genomes gff dataframes
    names_columns_gff = ["contig","source","type", "left_coordinate", "right_coordinate", "dot", "strand", "0" ,"description"] 
    d_gff = {}
    for g in list_genomes:
        d_gff[g] = pd.read_csv(f"{protein_gff_folder}/gff3/{g}.gff", comment = "#", sep = "\t", names = names_columns_gff)
        d_gff[g]["description"] = d_gff[g]["description"].apply(lambda x: x.split(';')[0].split('=')[-1])

    # First we create a nested dict with k = TA, v = {genome: [list of list_genes_same_corespot]}
    d_corespot_genes2test = {}
    df_core_features = pd.read_csv(core_features_file, sep = "\t", header = 0, comment ="#", index_col = 0)
    for TA, list_genome_to_search_in in d_genome_with_noTA_same_core_spot.items():
        d_corespot_genes2test[TA] = {}
        for g in list_genome_to_search_in:
            tmp_core1_row = df_core_features[(df_core_features["core_family"] == int(d_TA_core_ref[TA][0].rsplit(".",1)[0])) & (df_core_features["genome_name"] == g)]
            tmp_core2_row = df_core_features[(df_core_features["core_family"] == int(d_TA_core_ref[TA][1].rsplit(".",1)[0])) & (df_core_features["genome_name"] == g)]
            # Easiest case, both core are on the same contig, in this case we keep the genes
            if tmp_core1_row["contig"].values[0] == tmp_core2_row["contig"].values[0]:
                tmp_index1 = d_gff[g][d_gff[g]["description"] == tmp_core1_row.index[0]].index[0]
                tmp_index2 = d_gff[g][d_gff[g]["description"] == tmp_core2_row.index[0]].index[0]
                d_corespot_genes2test[TA][g] = [d_gff[g].iloc[min(tmp_index1,tmp_index2)+1:max(tmp_index1,tmp_index2)]["description"].tolist()]
            
            # In case, we have the two core genes on different contigs in this genome, we first try to determine if there are others core genes within the same contig 
            # which would help us to determine the probable contig orientation
            elif tmp_core1_row["contig"].values[0] != tmp_core2_row["contig"].values[0] :
                #then it exists 4 possible combinaison of the core spot reconstitution
                list_of_list_genes2test = []
                tmp_list_1 = []
                tmp_list_2 = []
                tmp_list_3 = []
                tmp_list_4 = []

                tmp_index1 = d_gff[g][d_gff[g]["description"] == tmp_core1_row.index[0]].index[0]
                tmp_index2 = d_gff[g][d_gff[g]["description"] == tmp_core2_row.index[0]].index[0]
                #there are others core genes within these contigs
                if len(df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]]) >=2:
                    #and the core gene is either the first one or the last one on this contig (then we take the genes on the extremity of the contig)
                    if df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]].sort_values(by = "left_coordinate").index[0] == tmp_core1_row.index[0] :
                        tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core1_row["left_coordinate"].values[0]))]["description"].tolist()

                    elif df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]].sort_values(by = "left_coordinate").index[-1] == tmp_core1_row.index[0] :
                        tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core1_row["left_coordinate"].values[0]))]["description"].tolist()

                    # in case the core gene is not at the extremity of the contig (maybe due to wrong assembly/ large rearragement?) 
                    # we keep the genes up to the next core genes (2 combinaisons then)
                    elif (df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]].sort_values(by = "left_coordinate").index[-1] != tmp_core1_row.index[0]) and df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]].sort_values(by = "left_coordinate").index[0] != tmp_core1_row.index[0]:
                        # need to get the index of the two closest core from our core reference gene
                        #it is a bit messy, will try later to change that
                        df_core_this_contig = df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]].sort_values(by = "left_coordinate").reset_index(drop = False)
                        new_neighbouring_core_genes = [df_core_this_contig.iloc[df_core_this_contig[df_core_this_contig["index"] == tmp_core1_row.index[0]].index[0]-1],
                                                    df_core_this_contig.iloc[df_core_this_contig[df_core_this_contig["index"] == tmp_core1_row.index[0]].index[0]+1]]

                        tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core1_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()
                        tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core1_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()


                # if there are no others core genes on the same contig, we do create both combinaisons
                elif len(df_core_features[df_core_features["contig"] == tmp_core1_row["contig"].values[0]]) <2:
                    tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core1_row["left_coordinate"].values[0]))]["description"].tolist()    
                    tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core1_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core1_row["left_coordinate"].values[0]))]["description"].tolist()

                # now looking for the second reference core gene
                # Note: we need to be careful now that there are already two possible existing list of genes
                if len(df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]]) >=2:
                    if df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]].sort_values(by = "left_coordinate").index[0] == tmp_core2_row.index[0] :
                        if tmp_list_1 != []:
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        if tmp_list_2 != []:
                            tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        if tmp_list_1 == [] and tmp_list_2 == []:
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()


                    elif df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]].sort_values(by = "left_coordinate").index[-1] == tmp_core2_row.index[0] :
                        if tmp_list_1 != []:
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        if tmp_list_2 != []:
                            tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        if tmp_list_1 == [] and tmp_list_2 == []:
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()


                    #finally if the second core gene is also in the middle of new core genes we have 4 combinaisons
                    elif (df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]].sort_values(by = "left_coordinate").index[-1] != tmp_core2_row.index[0]) and df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]].sort_values(by = "left_coordinate").index[0] != tmp_core2_row.index[0]:
                        df_core_this_contig = df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]].sort_values(by = "left_coordinate").reset_index(drop = False)
                        new_neighbouring_core_genes = [df_core_this_contig.iloc[df_core_this_contig[df_core_this_contig["index"] == tmp_core2_row.index[0]].index[0]-1],
                                                    df_core_this_contig.iloc[df_core_this_contig[df_core_this_contig["index"] == tmp_core2_row.index[0]].index[0]+1]]

                        if tmp_list_1 != []:
                            if d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist() != [] and d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist() != []:
                                tmp_list_3 = tmp_list_1 + d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()
                                tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()
                                
                            elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist() != [] :
                                tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()

                            elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist() != [] :
                                tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()

                        if tmp_list_2 != []:
                            if d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist() != [] and d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist() != []:
                                tmp_list_4 = tmp_list_2 + d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()
                                tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()
                                
                            elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist() != [] :
                                tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()

                            elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist() != [] :
                                tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()

                        if tmp_list_1 == [] and tmp_list_2 == [] :
                            tmp_list_1 = d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] > int(new_neighbouring_core_genes[0]["left_coordinate"]))]["description"].tolist()
                            tmp_list_2 = d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0])) & (d_gff[g]["left_coordinate"] < int(new_neighbouring_core_genes[1]["left_coordinate"]))]["description"].tolist()

                #if the second core gene is the only one within its contig
                elif len(df_core_features[df_core_features["contig"] == tmp_core2_row["contig"].values[0]]) <2:
                    if tmp_list_1 != []:
                        if d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != [] and d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != [] :
                            tmp_list_3 = tmp_list_1 + d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        
                        elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != []:
                            tmp_list_1 =+ d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        
                        elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != []: 
                            tmp_list_1 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()

                    elif tmp_list_2 != [] :
                        if d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != [] and d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != [] :
                            tmp_list_4 = tmp_list_2 + d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                            tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        
                        elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != []:
                            tmp_list_2 =+ d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        
                        elif d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist() != []: 
                            tmp_list_2 += d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()

                    elif tmp_list_1 == [] and tmp_list_2 == [] :
                        tmp_list_1 = d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] < int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()
                        tmp_list_2 = d_gff[g][(d_gff[g]["contig"] == tmp_core2_row["contig"].values[0]) & (d_gff[g]["type"] == "CDS") & (d_gff[g]["left_coordinate"] > int(tmp_core2_row["left_coordinate"].values[0]))]["description"].tolist()

                for number_tmp_list in [tmp_list_1, tmp_list_2, tmp_list_3, tmp_list_4]:
                    if number_tmp_list != []:
                        list_of_list_genes2test.append(number_tmp_list)

                d_corespot_genes2test[TA][g] = list_of_list_genes2test

    list_columns = ["qseqid", "qlen", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    str_columns_blastp = "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    # now doing the blastp of the reference neighbouring genes with all genes of all possibles combinaisons of every genomes for each TA
    ref_neigbours_seqio = SeqIO.index(f"{tmp_folder}/ref_genome_neighbours_db.faa", "fasta")
    #now create a query_tmp file with the sequences of the ref neighbouring genes
    for TA in df_all_neigbouring_genes["TA_homolog_of"].drop_duplicates().tolist():
        list_ref_neigbours_this_TA = df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == reference_genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)].drop_duplicates(subset = "TA_homolog_of")["neighbours_genes"].values[0]
        list_ref_neigbours_this_TA = [x for x in list_ref_neigbours_this_TA.split(",") if x]
        with open(f"{tmp_folder}/TA_neigbouring_genes_tmp_blastp_query.faa", "w") as tmp_query:
            for i in list_ref_neigbours_this_TA:
                tmp_query.write(f">{ref_neigbours_seqio[i].id}\n{ref_neigbours_seqio[i].seq}\n")

        # creating a tmp list file with all genes ids to include in the blastp search
        with open(f"{tmp_folder}/list_allgenes_samecore_noTA.lst", "w") as f :
            for list_of_list2_combine in d_corespot_genes2test[TA].values():
                for l in list_of_list2_combine:
                    for gene in l:
                        f.write(f"{gene}\n")
    
        # blastp
        blastp_cline = NcbiblastpCommandline(query = f"{tmp_folder}/TA_neigbouring_genes_tmp_blastp_query.faa", db = f"{tmp_folder}/all_proteome_db", evalue = 0.1,
                                            outfmt = str_columns_blastp , seqidlist = f"{tmp_folder}/list_allgenes_samecore_noTA.lst",
                                            out= f"{outdir}/2-neighbouring_genes/blastp_raw_data/{TA}_blastp_neighbours_within_corespot_noTA_all_combinaisons.csv")
        stdout, stderr = blastp_cline()

    return d_corespot_genes2test, d_genome_with_TA_same_core_spot


def analyze_blastp_neighbouring_genes_step(df_all_neigbouring_genes, d_corespot_testedgenes_noTA, reference_genome, list_genomes, outdir):

    list_genomes = [x.rsplit(".",1)[0] for x in list_genomes if x]
    n_index_final_df = 0
    d_res = {}
    list_columns = ["qseqid", "qlen", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df_TA_tblast = pd.read_csv(f"{outdir}/1-tblastn/full_TA_tblastn_with_core.csv", sep = "\t", index_col= 0)

    for TA, d_genes_tested_in_each_genome in d_corespot_testedgenes_noTA.items():
        df_neighbours_near_TA = pd.read_csv(f"{outdir}/2-neighbouring_genes/blastp_raw_data/{TA}_blastp_neighbours_near_TA.csv", sep = "\t", names= list_columns)
        df_neighbours_noTA = pd.read_csv(f"{outdir}/2-neighbouring_genes/blastp_raw_data/{TA}_blastp_neighbours_within_corespot_noTA_all_combinaisons.csv", sep = "\t", names= list_columns)

        for genome in list_genomes:
            #first we reset every variable we store in our final dict
            is_TA_present = "-"
            is_same_corespot = "-"
            TA_multiple_copies = "-"
            TA_contig = "-"
            corefamily_left = "-"
            corefamily_right = "-"
            tox_pid = "-"
            tox_pcov ="-"
            tox_left_c ="-"
            tox_right_c ="-"
            antitox_pid = "-"
            antitox_pcov = "-"
            antitox_left_c = "-"
            antitox_right_c = "-"
            list_neighbours_TA_this_genome = "-"
            list_original_corespot_conserved_genes = "-"
            tuple_number_original_corespot_conserved_compared_to_ref = ("-","-")
            TA_different_location = "-"
            tuple_number_different_location_conserved_compared_to_ref = ("-","-")

            if genome == reference_genome:
                is_reference_genome = "Yes"
                ref_corefamily_left = df_all_neigbouring_genes[df_all_neigbouring_genes["TA_homolog_of"] == TA]["left_core_family"].tolist()[0]
                ref_corefamily_right = df_all_neigbouring_genes[df_all_neigbouring_genes["TA_homolog_of"] == TA]["right_core_family"].tolist()[0]
                corefamily_left = ref_corefamily_left
                corefamily_right = ref_corefamily_right
                number_neighbours_genes_ref = df_all_neigbouring_genes[df_all_neigbouring_genes["TA_homolog_of"] == TA]["n_neighbours"].tolist()[0]
                list_ref_neighbours = [x for x in df_all_neigbouring_genes[df_all_neigbouring_genes["TA_homolog_of"] == TA]["neighbours_genes"][0].split(",") if x]

            else :
                is_reference_genome = "No"

            if df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])].empty == False :
                is_TA_present = "Yes"
                # if there are only one copy of the TA system on the studied genome
                if len(df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]) == 1:
                    TA_multiple_copies = "No"

                    corefamily_left = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["left_core_family"].values[0]
                    corefamily_right = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["right_core_family"].values[0]
                    is_same_corespot = check_corespot_with_ref(corefamily_left, corefamily_right, ref_corefamily_left,  ref_corefamily_right) 
                    
                    if is_same_corespot == "Yes" or is_same_corespot == "?":
                        TA_contig = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["sseqid_TA"].values[0]
                        tox_pid = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pident"].values[0]
                        tox_pcov = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pcov"].values[0]
                        tox_left_c = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_left_coordinate"].values[0]
                        tox_right_c = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_right_coordinate"].values[0]
                        antitox_pid = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pident"].values[0]
                        antitox_pcov = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pcov"].values[0]
                        antitox_left_c = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_left_coordinate"].values[0]
                        antitox_right_c = df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_right_coordinate"].values[0]
                        list_neighbours_TA_this_genome = [x for x in df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]["neighbours_genes"].values[0].split(",") if x]
                        list_original_corespot_conserved_genes = check_neighbouring_genes_conservation(list_neighbours_TA_this_genome, list_ref_neighbours, df_neighbours_near_TA)
                        if genome == reference_genome :
                            tuple_number_original_corespot_conserved_compared_to_ref = (len(list_ref_neighbours), len(list_ref_neighbours))
                        else :
                            tuple_number_original_corespot_conserved_compared_to_ref = (len(list_original_corespot_conserved_genes), len(list_ref_neighbours))

                    elif is_same_corespot == "No":
                        TA_different_location = [df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["sseqid_TA"].values[0],
                                                            min(df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_left_coordinate"].values[0],df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_left_coordinate"].values[0]),
                                                            max( df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_right_coordinate"].values[0], df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_right_coordinate"].values[0]),
                                                            df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pident"].values[0],
                                                            df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pcov"].values[0],
                                                            df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pident"].values[0],
                                                            df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pcov"].values[0]
                                                            ]
                        
                        list_TA_different_location_neighbours = [x for x in df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]["neighbours_genes"].values[0].split(",") if x]
                        list_TA_different_location_neighbours_conserved = check_neighbouring_genes_conservation(list_TA_different_location_neighbours, list_ref_neighbours, df_neighbours_near_TA)
                        tuple_number_different_location_conserved_compared_to_ref = (len(list_TA_different_location_neighbours_conserved), len(list_ref_neighbours))


                        list_original_corespot_conserved_genes = best_combinaison_of_genes_conserved(TA, genome, d_corespot_testedgenes_noTA, df_neighbours_noTA, list_ref_neighbours)
                        tuple_number_original_corespot_conserved_compared_to_ref = (len(list_original_corespot_conserved_genes), len(list_ref_neighbours))
        
                elif len(df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]) > 1:
                    TA_multiple_copies = "Yes"
                    n_break = 2
                    is_same_corespot = "No"

                    for row in df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])].iterrows():
                        n_break -=1
                        if n_break == 0:
                            break

                        tmp_corefamily_left = row[1]["left_core_family"]
                        tmp_corefamily_right = row[1]["right_core_family"]
                        tmp_is_same_corespot = check_corespot_with_ref(corefamily_left, corefamily_right, tmp_corefamily_left,  tmp_corefamily_right) 

                        if tmp_is_same_corespot == "Yes" or is_same_corespot == "?":
                            is_same_corespot = "Yes"
                            corefamily_left = tmp_corefamily_left
                            corefamily_right = tmp_corefamily_right
                            TA_contig = row[1]["sseqid_TA"]
                            tox_pid = row[1]["tox_pident"]
                            tox_pcov = row[1]["tox_pcov"]
                            tox_left_c = row[1]["tox_left_coordinate"]
                            tox_right_c = row[1]["tox_right_coordinate"]
                            antitox_pid = row[1]["antitox_pident"]
                            antitox_pcov = row[1]["antitox_pcov"]
                            antitox_left_c = row[1]["antitox_left_coordinate"]
                            antitox_right_c = row[1]["antitox_right_coordinate"]
                            list_neighbours_TA_this_genome = [x for x in df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]["neighbours_genes"].values[0].split(",") if x]
                            list_original_corespot_conserved_genes = check_neighbouring_genes_conservation(list_neighbours_TA_this_genome, list_ref_neighbours, df_neighbours_near_TA)
                            if genome == reference_genome:
                                tuple_number_original_corespot_conserved_compared_to_ref = (len(list_original_corespot_conserved_genes), len(list_original_corespot_conserved_genes))
                            else :
                                tuple_number_original_corespot_conserved_compared_to_ref = (len(list_original_corespot_conserved_genes), len(list_ref_neighbours))

                        elif tmp_is_same_corespot == "No":
                            TA_different_location = [df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["sseqid_TA"].values[0],
                                                                min(df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_left_coordinate"].values[0],df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_left_coordinate"].values[0]),
                                                                max( df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_right_coordinate"].values[0], df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_right_coordinate"].values[0]),
                                                                df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pident"].values[0],
                                                                df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["tox_pcov"].values[0],
                                                                df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pident"].values[0],
                                                                df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])]["antitox_pcov"].values[0]
                                                                ]
                            
                            list_TA_different_location_neighbours = [x for x in df_all_neigbouring_genes[(df_all_neigbouring_genes["genome_name"] == genome) & (df_all_neigbouring_genes["TA_homolog_of"] == TA)]["neighbours_genes"].values[0].split(",") if x]
                            list_TA_different_location_neighbours_conserved = check_neighbouring_genes_conservation(list_TA_different_location_neighbours, list_ref_neighbours, df_neighbours_near_TA)
                            tuple_number_different_location_conserved_compared_to_ref = (len(list_TA_different_location_neighbours_conserved), len(list_ref_neighbours))


            elif df_TA_tblast[(df_TA_tblast["genome_name"] == genome) & (df_TA_tblast["ref_toxin"] == TA.split("-")[0])].empty == True :
                is_TA_present = "No"
                list_original_corespot_conserved_genes = best_combinaison_of_genes_conserved(TA, genome, d_corespot_testedgenes_noTA, df_neighbours_noTA, list_ref_neighbours)
                tuple_number_original_corespot_conserved_compared_to_ref = (len(list_original_corespot_conserved_genes), len(list_ref_neighbours))

        
            #now adding the results to the dictionnary

            if tuple_number_different_location_conserved_compared_to_ref == ("-","-"):
                n_conserved_diff_core = "-"
                percentage_conserved_diff_core = "-"
            else :
                n_conserved_diff_core = f"{tuple_number_different_location_conserved_compared_to_ref[0]}/{tuple_number_different_location_conserved_compared_to_ref[1]}"
                percentage_conserved_diff_core = int(tuple_number_different_location_conserved_compared_to_ref[0])/int(tuple_number_different_location_conserved_compared_to_ref[1])*100
            
            if tuple_number_original_corespot_conserved_compared_to_ref == ("-","-"):
                n_conserved_same_core = "-"
                percentage_conserved_same_core = "-"
            else :
                n_conserved_same_core = f"{tuple_number_original_corespot_conserved_compared_to_ref[0]}/{tuple_number_original_corespot_conserved_compared_to_ref[1]}"
                percentage_conserved_same_core = int(tuple_number_original_corespot_conserved_compared_to_ref[0])/int(tuple_number_original_corespot_conserved_compared_to_ref[1])*100

            d_res[n_index_final_df] =   {"Ref_TA": TA,
                                        "Genome": genome,
                                        "Is_TA_homolog": is_TA_present,
                                        "Same_core_spot": is_same_corespot,
                                        "Multiples_copies": TA_multiple_copies,
                                        "TA_contig": TA_contig,
                                        "Core_family_left": corefamily_left,
                                        "Core_family_right": corefamily_right,
                                        "tox_pid": tox_pid,
                                        "tox_pcov" : tox_pcov,
                                        "tox_left_coordinate": tox_left_c,
                                        "tox_right_coordinate": tox_right_c,
                                        "antitox_pid": antitox_pid,
                                        "antitox_pcov": antitox_pcov,
                                        "antitox_left_coordinate": antitox_left_c,
                                        "antitox_right_coordinate": antitox_right_c,
                                        "list_neighbours_TA": list_neighbours_TA_this_genome,
                                        "list_conserved_compared2ref": list_original_corespot_conserved_genes,
                                        "n_conserved_genes/ref_genes": n_conserved_same_core,
                                        "%_conserved_genes/ref_genes": percentage_conserved_same_core,
                                        "TA_other_corespot": TA_different_location,
                                        "n_conserved_diff_core": n_conserved_diff_core,
                                        "%_conserved_diff_core": percentage_conserved_diff_core
                                        }

            n_index_final_df += 1

    df_res = pd.DataFrame.from_dict(d_res, orient = "index")
    return df_res

def check_corespot_with_ref(core1, core2, ref_core1, ref_core2):
    #small function which return "Yes" if the core genes in the new genomes match with the reference corespot genes family and "No" if not the case
    #Note: in case the TA in the new genome is within a contig with only one core gene, the other one would be "None"
    #In case the unique core here match with one of the ref core gene it will return "Yes"
    #However if the two coregenes are known ( and so different of None), both needs to match with the two ref core genes to return "Yes"
    #if the TA is within a contig without any core genes => it will return "?"

    if core1 != "None" and core2 != "None":
        if core1 in [ref_core1, ref_core2] and core2 in [ref_core1, ref_core2]:
            return "Yes"
        else :
            return "No"

    elif core1 == "None" and core2 == "None":
        return "?"
    
    elif core1 != "None":
        if core1 in [ref_core1, ref_core2]:
            return "Yes"
        else:
            return "No"
    
    elif core2 != "None":
        if core2 in [ref_core1, ref_core2]:
            return "Yes"
        else:
            return "No"


def check_neighbouring_genes_conservation(list_genes_new_genome, list_genes_ref, df_blastp_results):
    #return a list of genes which matched through blastp (threshold 80% pid, 80% pcov)
    df_hits_neighbouring_genes = df_blastp_results[(df_blastp_results["pident"] > 80) 
                                                    & (df_blastp_results["sseqid"].isin(list_genes_new_genome)) & (df_blastp_results["qseqid"].isin(list_genes_ref))
                                                    & (df_blastp_results["length"]/df_blastp_results["qlen"]*100 > 80)
                                                    ].sort_values(by ="evalue")
    res_list = []
    for gene in list_genes_ref:
        if df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene].empty == False :
            if len(df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene]) == 1 and df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene]["sseqid"].values[0] not in res_list:
               res_list.append(df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene]["sseqid"].values[0])
            elif len(df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene]) > 1 :
                for gene_matching_with_ref in df_hits_neighbouring_genes[df_hits_neighbouring_genes["qseqid"] == gene]["sseq"].tolist():
                    if gene_matching_with_ref not in res_list:
                        res_list.append(gene_matching_with_ref)
                        break
    
    return res_list

def best_combinaison_of_genes_conserved(TA_name, genome_name, d_corespot_testedgenes_noTA, df_blastp_results, list_genes_ref):
    # because there may be up to 4 combinaisons to test, we determine which one of these combinaison give the best results
    best_combinaison = []
    max_len = 1000000 #just so we can select the minimal combinaison of genes with the highest number of conserved genes
    for gene_combinaison_list in d_corespot_testedgenes_noTA[TA_name][genome_name]:
        tmp_res = check_neighbouring_genes_conservation(gene_combinaison_list, list_genes_ref, df_blastp_results)
        if len(tmp_res) > len(best_combinaison):
            best_combinaison = tmp_res
            max_len = len(gene_combinaison_list)
        elif len(tmp_res) == len(best_combinaison) and len(tmp_res) != 0 and len(gene_combinaison_list) < max_len:
            best_combinaison = tmp_res
            max_len = len(gene_combinaison_list)

    return best_combinaison



####################################################################################################################################

class MultiIndexDict:
    def __init__(self, *indexes):
        self._indexes = indexes
    def __getitem__(self, key):
        for idx in self._indexes:
            try:
                return idx[key]
            except KeyError:
                pass
        raise KeyError("{0} not found".format(key))


if __name__ == "__main__" :
        main()