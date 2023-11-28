from neuromaps.datasets import fetch_atlas
from neuromaps.datasets import available_annotations
from neuromaps.datasets import fetch_annotation
import nibabel as nib
import os
from os.path import expanduser
# set home path
home = expanduser("~")


# get fsLR
fslr = fetch_atlas('fsLR', '32k')

# get all available annotations
all_annotations=available_annotations()

# get first element of every sublist 
# https://www.geeksforgeeks.org/python-get-first-element-of-each-sublist/
def Extract(lst):
    return [item[0] for item in lst]
    
all_annot=Extract(all_annotations)

# download annotations, and continue if error for an individual atlas
annotations=[]
for annot in all_annot:
	try:
		print(annot)
		annotation = fetch_annotation(source=annot, return_single=False)
		annotations.append(annotation)
	except:
		pass
		
# get list of downloaded annotations
annot_dirs = os.listdir(home+"/neuromaps-data/annotations")

# manually remove some maps that didn't fully download
annot_use = annot_dirs
annot_use.sort()
# raichle cmrglc didn't download
annot_use.remove("raichle")

# re-fetch annotations to get dictionary with only downloaded maps
for i in range(0, len(annot_use)):
	print(annot_use[i])
	annot_have = fetch_annotation(source=annot_use[i], return_single=False)	