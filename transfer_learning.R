# Transfer learning to generate image features by subward (after generating patches for each subward)

library(keras)
library(tensorflow)
library(EBImage)
library(stringr)

patch_dir <- "subward_patches"

transfer_model <- application_inception_resnet_v2
patch_dims <- c(200, 200, 3)

conv_base <- transfer_model( 
  weights="imagenet",
  include_top=FALSE,
  input_shape=patch_dims
)

patch_file_list <- list.files(patch_dir)
features_matrix <- numeric(0)
batch_size <- 100 # feed patches into transfer model 100 at a time
nbatches <- ceiling(length(patch_file_list) / batch_size)

for ( i in 1:nbatches ){
  
  # List of patches to include in batch
  if ( i == nbatches ){
    batch_size <- length(patch_file_list) - (i - 1) * batch_size # as there may be fewer than batch_size patches left
  }
  batch_patches <- patch_file_list[((i - 1) * batch_size + 1):(i * batch_size)]
  
  # Read in files and put them in an array
  patch_array <- array(0, dim=c(batch_size, patch_dims))
  
  k <- 0
  for ( filename in batch_patches ){
    k <- k + 1
    patch <- readRDS(paste0(patch_dir, "/", filename))
    patch_array[k,,, ] <- patch
  }
  
  # Get features
  features_array <- conv_base %>% predict(patch_array)
  
  # Maxpooling
  features_array_maxpooled <- array(0, dim=dim(features_array)[c(1, 4)])
  for ( dim1 in 1:dim(features_array)[1] ){
    for ( dim2 in 1:dim(features_array)[4] ){
      features_array_maxpooled[dim1, dim2] <- max(features_array[dim1,,,dim2])
    }
  }
  features_matrix <- rbind(features_matrix, features_array_maxpooled)

  print(i)
}
rownames(features_matrix) <- str_replace_all(patch_file_list, ".rds", "")

#################################################################################################################################
# Get features for each subward

patch_features <- readRDS("Results new - April 2021/patch_features_high_dim_irv2_maxpooled_new.rds")

subward_indices <- 12363:12814
subward_features <- subward_list <- numeric(0)

for ( subward_index in subward_indices ){

  patch_features_subward <- patch_features[grepl(subward_index, rownames(patch_features)), ]

  if ( dim(patch_features_subward)[1] > 0 ){

    # Take means of features for subward patches
    subward_features <- rbind(subward_features, colMeans(patch_features_subward))
    subward_list <- c(subward_list, subward_index)

  }

}

rownames(subward_features) <- subward_list - 12362 # convert to subward_id's used in dar.shapefile
write.csv(subward_features, "subward_image_features.csv")
