# 22q_spatial_transcriptomic
Analyses relating spatial patterns of gene expression from Allen Human Brain Atlas (AHBA) to fMRI markers of 22q11.2 Deletion Syndrome (22qDel)

## Overview
Mouse model research suggests that deficits in specific neuronal cell types may underlie altered brain activity and behavior in humans with 22qDel. Specifically, multiple 22qDel mouse models converge on disruptions to parvalbumin expressing (PV+) inhibitory interneurons. There is also evidence from the 22qDel mouse model for disrupted long-range connectivity of excitatory pyramidal neurons (PN). In humans with 22qDel, PV+ and PN dysfunction may contribute to disrupted brain activity observable with fMRI, but direct invasive confirmation of this hypothesis is not feasible. Spatial transcriptomics provides a non-invasive means to relate human neuroimaging findings to typical patterns of brain gene expression. We can test the hypothesis that fMRI disruptions in 22qDel index regions normally enriched for specific neuronal marker genes (e.g. PVALB for PV+ interneurons). 

We will compute several measures from resting-state fMRI across a set of 718 cortical and subcortical parcels from the [CAB-NP atlas](https://github.com/ColeLab/ColeAnticevicNetPartition) -- specifically global brain connectivity (GBC), local connectivity/regional homogeneity (ReHo), and temporal variability (RSFA). Repository for first-level fMRI analyses can be found [here](https://github.com/cschleifer/22q_hoffman). Scripts in the current repository extract AHBA gene expression data from the same CAB-NP parcels using the [abagen toolbox](https://abagen.readthedocs.io/en/stable/usage.html), read these data along with parcellated fMRI measures, generate brain plots, and test spatial relationships between gene expression and fMRI group differences.  

## Dependencies
* Requires wb_command for ciftiTools functions that read/plot MRI data. 
  * Download: https://www.humanconnectome.org/software/get-connectome-workbench
  * Script expects the workbench directory to be `/Applications/workbench/bin_macosx64` (either download to this location, or edit the path in the script)
* To read fMRI results from the hoffman2 server, the script expects that the server `hoffman2.idre.ucla.edu:/u/project/cbearden/data` is mounted to your local machine at the path `~/Desktop/hoffman_mount` using an application such as SSHFS (mac download: https://osxfuse.github.io/)
  * Requires first-level MRI results to be already computed on server (see https://github.com/cschleifer/22q_hoffman)
  
## Analysis Steps
1. Generate parcellated BOLD metrics (e.g. RSFA) on hoffman
    * Scripts in this repo are intended to be run locally but require inputs files that are best generated on a HPC cluster (e.g. hoffman2) due to storage and computational constraints. For code to be run on the cluster prior to these scripts, see: https://github.com/cschleifer/22q_hoffman 


2. Convert CAB-NP parcellation in CIFTI space to fsaverage5 surface and MNI volume space for input into abagen
    * [prep_cabnp_abagen.sh](prep_cabnp_abagen.sh) uses wb_command to convert CIFTI dlabel file 
    * Outputs are in [CAB_NP_CONVERTED](CAB-NP/CAB_NP_converted) directory


3. Use abagen to parcellate AHBA data based on CAB-NP atlas
    * [abagen_ahba_cabnp.py](abagen_ahba_cabnp.py) uses abagen.get_expression_data to extract gene expression for each parcel
    * Outputs are [CAB-NP_subcort_abagen_expression.csv](CAB-NP_surface_abagen_expression.csv) and [CAB-NP_subcort_abagen_expression.csv](CAB-NP_subcort_abagen_expression.csv)
      * Surface and volume are computed separately -- surface using fsaverage5 space, volume using MNI. Expression is normalized separately within brain structures (cortex, cerebellum, subcort/brainstem) because computing volume and surface separately complicates normalizing across the whole brain 


4. Test relationship between PVALB/SST expression gradient and fMRI group difference (22qDel - healthy controls)
    * [22q_spatial_transcriptomic.R](22q_spatial_transcriptomic.R) reads parcellated AHBA data into R, along with subject level MRI measures in csv format (see Step 1)
    * Requires wb_command for ciftiTools functions that read/plot MRI data. 
      * Download: https://www.humanconnectome.org/software/get-connectome-workbench
      * Script expects the workbench directory to be `/Applications/workbench/bin_macosx64` (either download to this location, or edit the path in the script)
      * To read fMRI results from the hoffman2 server, the script expects that the hoffman2 server (`hoffman2.idre.ucla.edu:/u/project/cbearden/data`) is mounted to your local machine (e.g. using [SSHFS](https://osxfuse.github.io/)) at the path `~/Desktop/hoffman_mount`. Either mount at this path or adjust the `hoffman` variable in the script
 
