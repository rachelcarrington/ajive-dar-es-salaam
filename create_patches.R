# Generate random patches for each subward

library(EBImage)
library(stringr)

# Set directories
subward_image_dir <- "subward_images"
patch_dir <- "subward_patches"
dir.create(patch_dir)

# Get image file list
subward_list <- list.files(subward_image_dir)

# Set number of patches per subward and patch dimensions
num_patches_to_generate <- 10
patch_dims <- c(200, 200)

for ( i in 1:length(subward_list) ){
  
  # Read in subward image
  filename <- subward_list[i]
  dat <- readImage(paste0(subward_image_dir, "/", filename))
  
  # Create patch filename
  patch_filename_root <- str_replace(filename, ".tif", "")
  patch_filename_root <- str_replace(patch_filename_root, "subward_", "patch_")
  patch_filename_root <- paste0(patch_dir, "/", patch_filename_root)
  
  # Range of possible starting points (top left corners) for patches
  x_range <- 1:(dim(dat)[1] - (patch_dims[1] - 1))
  y_range <- 1:(dim(dat)[2] - (patch_dims[2] - 1))
  
  # Generate patches for subward
  npatches <- 0
  
  while( npatches < num_patches_to_generate ){
    
    # Generate patch
    x <- sample(x_range, 1)
    y <- sample(y_range, 1)
    patch <- dat[x:(x + 199), y:(y + 199), ]
    
    # If less than 1% background, save patch and increase npatches by 1
    if ( sum(patch == 0) <= 1200 ){
      npatches <- npatches + 1
      patch_filename <- paste0(patch_filename_root, "_", npatches, ".rds")
      saveRDS(patch, patch_filename)
    }
    
  }
  
  rm("dat") # to free up memory  
}
