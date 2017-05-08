# Brain Plots
#### Currently only supports openfmri group-level analysis

#### How to Use
* Nibabel master (after merging nibabel/#249) **Contains CIFTI2 support!**
* Clone repo (including 32k Atlas)
* [Recommended] Use commented command at the top to run windowless

##### Arguments
  * -d (__required__): path to group results (up to tasks)
  * -t (__required__): tasks separated by spaces
  * -o: path to output directory (default is current working directory)
  * -a: path to `32k_ConteAtlas_v2` folder (default is current working directory)
  * -i: flag to make uninflated brain output (default is inflated)
  * -s: flag to make whole brain output (default is split)
  * -ss: flag to remove quad split
  * -th: set threshold (must be a float) (default is 2.3)
  * -dt: set display threshold (int) (default is 6)
 
##### Sample Use
    om_py_plot_surfs.py -d /om/project/voice/processedData/groupAnalysis/may2016_model101/one_sample/model101 -t task001 task002 -ss -th 3.0 -dt 5
 


### TO-DO
- [ ] Add `mask` input to show brain coverage
- [ ] Mutable image size
- [ ] Handle flexible directory structure - Subject Level QA, Group Analysis, Etc
