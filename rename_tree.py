import sys


#small script which rename the terminals nodes of a phylogenetical tree given a namefile as follow
# one line : {original_name}:{new_name}

#commandline example : python script/rename_tree.py panacota_redo_photo_morga/22-12-12_photo_mm.treefile panacota_redo_photo_morga/redo_name_dict.txt


def main():

    tree2change = sys.argv[1]
    namefile = sys.argv[2]

    d_name = {}
    with open(namefile, "r") as f:
        for line in f:
            d_name[line.split(":",1)[0].strip()] = line.split(":",1)[1].strip()
    

    with open(tree2change,"r") as f:
        str_new_tree = f.read()
    

    for k,v in d_name.items():
        str_new_tree = str_new_tree.replace(k,v)
    
    with open(tree2change.rsplit(".",1)[0]+"_renamed.nwk", "w") as f:
        f.write(str_new_tree)


if __name__ == "__main__":
    main()