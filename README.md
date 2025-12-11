# Urban mapping in Dar es Salaam using AJIVE

This repository contains code and data to accompany the paper [**Urban mapping in Dar es Salaam using Angle-Based Joint and Individual Variation Explained**](https://academic.oup.com/jrsssc/advance-article/doi/10.1093/jrsssc/qlaf043/8236692).

The AJIVE algorithm we used was developed by [(Feng et al., 2018)](https://www.sciencedirect.com/science/article/pii/S0047259X1730204X).

The satellite image data used in the paper can be downloaded from [MAXAR](https://www.maxar.com/open-data/covid19).
The subward-level features used in our analysis are available here in the file subward_features.csv

## Functions

* ajive.R : implement AJIVE algorithm. If the joint rank is 0, implements PCA on data blocks individually.
	(need to fix it to estimate initial ranks if not supplied)
* create_patches.R : code used to create random patches from images
* estimate_joint_rank.R : estimate the joint rank of 2 or more data matrices
* jive.R : implementation of the JIVE algorithm (Lock et al., 2013)[https://pmc.ncbi.nlm.nih.gov/articles/PMC3671601/]
* transfer_learning.R : implementing transfer learning on image data.
