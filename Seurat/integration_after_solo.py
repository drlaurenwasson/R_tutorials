#2021_09_17 Collect and integrate solo scores from files that aready have had solo run on them.
#Usage: python3 integration_after_solo.py <directory> <sample identifier folders>
import numpy as np
import scanpy as sc
import getopt, sys

argumentList = sys.argv
print(argumentList)
directory = sys.argv[1]
print("directory:" + directory)
identifiers = sys.argv[2:]
print("identifiers:" + str(identifiers))
identifiers = list(identifiers)
def solo_write_file_2(directory, identifier):
    softmax_scores = directory + "/" + identifier + "/softmax_scores.npy"
    is_doublet = directory + "/" + identifier + "/is_doublet.npy"
    h5ad= directory + "/" + identifier + "/" + identifier + "_pbmc.h5ad"
    print(softmax_scores)
    print(is_doublet)
    print(h5ad)
    solo1 = np.load(softmax_scores)
    solo2 = np.load(is_doublet)
    x=sc.read_h5ad(h5ad)
    
    x.obs['solo_score']=solo1
    x.obs['predicted_doublets_solo']=solo2
    x.obs.to_csv(identifier+"_pbmc_withsolo_METADATA.csv")

for identifier in identifiers:
    solo_write_file_2(directory, identifier)

filenames = []

for identifier in identifiers:
    filenames.append(identifier +"_pbmc_withsolo_METADATA.csv")
    print(filenames)

# Open file3 in write mode
with open('solo_all.txt', 'w') as outfile:
  
    # Iterate through list
    for names in filenames:
  
        # Open each file in read mode
        with open(names) as infile:
  
            # read the data from file1 and
            # file2 and write it in file3
            outfile.write(infile.read())
  
        # Add '\n' to enter data of file2
        # from next line
        outfile.write("\n")
