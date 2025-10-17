# Urban mapping in Dar es Salaam using AJIVE

This repository contains code and data to accompany the paper [**Urban mapping in Dar es Salaam using AJIVE**](https://arxiv.org/abs/2403.09014).

The AJIVE algorithm was developed by (Feng et al., 2018), [**Angle-Based Joint and Individual Variation Explained**](https://www.sciencedirect.com/science/article/pii/S0047259X1730204X).

The satellite image data used in the paper can be downloaded from MAXAR [here](*add link*).

## Functions

* ajive.R : implement AJIVE algorithm. If the joint rank is 0, implements PCA on data blocks individually.
	(need to fix it to estimate initial ranks if not supplied)
* create_patches.R : code used to create random patches from images
* 
