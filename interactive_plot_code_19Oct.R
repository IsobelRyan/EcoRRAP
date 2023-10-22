# -----------------------------------------------------------------------------
# LOAD LIBRARIES --------------------------------------------------------------
# -----------------------------------------------------------------------------

#install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")
#install.packages("dartR")
#install.packages("vcfR")

#library(gdsfmt)
#library(SNPRelate)
#library(vcfR)
#install.packages("adegenet")

library(adegenet)
library(dartR)
library(ggplot2)
library(tidyverse)
library(htmlwidgets)
library(plotly)
library(readxl)
library(magick) ## for resizing images
library(plotly) ## adding images

# -----------------------------------------------------------------------------
# LOAD DATA -------------------------------------------------------------------
# -----------------------------------------------------------------------------

# load in metadata
Stylo_clusters_metadata <- read_xlsx("Spis_taxa_metadata_13Oct_GRAPH.xlsx")

# load in pca data: zoe's code
Stylo_genlight <- gl.read.vcf("Stylo_sf095_noclones_ld.vcf")
Stylo_genlight$pop <- as.factor(Stylo_clusters_metadata$New_sample_name) 
Stylo_genlight$ind.names <- Stylo_clusters_metadata$New_sample_name


PCA_all <- glPca(Stylo_genlight, parallel= TRUE) 
PCA_all.df <- as.data.frame(PCA_all$scores)
PC=1:10
pve_all <- data.frame(PC, 100*PCA_all$eig[1:10]/sum(PCA_all$eig[1:10]))


# -----------------------------------------------------------------------------
# REMOVING AND RENAMING IMAGES IN FOLDER --------------------------------------
# -----------------------------------------------------------------------------

# specify directory where images are stored 
#directory <- "Spis_images"
#
## list the files in the directory with the respective pattern
#mid <- list.files(directory, pattern = "_MID[0-9]*.JPG", full.names = TRUE)
#macro <- list.files(directory, pattern = "_MACRO[0-9]*.JPG", full.names = TRUE)
#
## remove files with the respective pattern 
#file.remove(mid, macro)
#
#
## once the colony images are left, now want to remove _COLONY and any digits so that only the #sample_name is left 
#files <- list.files(directory, full.names = TRUE)
#
#for (file in files) {
#  new_name <- gsub("_COLONY[0-9]*", "", basename(file)) # also need to change CBHE to HER
#  new_path <- file.path(dirname(file), new_name)
#  file.rename(file, new_path)
#  cat("Renamed:", file, "to", new_path, "\n")
#}


# -----------------------------------------------------------------------------
# FOR RESIZING IMAGES ---------------------------------------------------------
# -----------------------------------------------------------------------------

#for (i in 1:nrow(PCA_all.df)) {
#  # Get the file path
#  file_path <- PCA_all.df$imagefile[i]
#  
#  # Load the image
#  img <- image_read(file_path)
#  
#  # Resize the image to 20%
#  resized_img <- image_resize(img, geometry = "20%")
#  
#  # Extract the file name without extension
#  file_name <- tools::file_path_sans_ext(basename(file_path))
#  
#  # Define the path to save the resized image
#  save_path <- paste0("Spis_images/resize/", file_name, ".JPG")
#  
#  # Save the resized image
#  image_write(resized_img, path = save_path)
#  
#  # Print message to indicate progress
#  cat("Processed:", file_path, "\n")
#}


# -----------------------------------------------------------------------------
# CREATE COLUMN FOR IMAGE FILE PATHS ------------------------------------------
# -----------------------------------------------------------------------------

# Assign Sample_name as row names
PCA_all.df$Sample_name <- rownames(PCA_all.df)

# Rename sample names to the original location names (this matches the New_sample_name column in metadata)
PCA_all.df$Sample_name <- gsub("CBHE", "HER", PCA_all.df$Sample_name)
PCA_all.df$Sample_name <- gsub("ONLI", "LIZ", PCA_all.df$Sample_name)
PCA_all.df$Sample_name <- gsub("TSAU", "AUK", PCA_all.df$Sample_name)
PCA_all.df$Sample_name <- gsub("ONMO", "MOO", PCA_all.df$Sample_name)


# Assign the folder where the images are stored
image_folder <- "Spis_images/resize"  

# Extract the file names for each image from the image folder 
image_filenames <- list.files(image_folder, pattern = "\\.JPG$", full.names = TRUE)

# Create the image URLs using the image filenames and the image folder
image_urls <- paste0(image_folder, "/", basename(image_filenames))

# Function to extract individual ID from the filenames
extract_id <- function(filename) {
  number <- sub("^[^0-9]*([0-9]+).*", "\\1", basename(filename))
  return(number)
}

# Use function to extract the ID
ID <- sapply(image_filenames, extract_id)

# Group the image filenames by their sample numbers
image_filenames_grouped <- split(image_filenames, ID)

# Create empty vector to assign image URLs
customdata <- character(length(PCA_all.df$Sample_name))

# Loop to assign the first image filename to each sample
for (i in seq_along(PCA_all.df$Sample_name)) {
  sample_number <- extract_sample_number(PCA_all.df$Sample_name[i])
  if (sample_number %in% names(image_filenames_grouped)) {
    customdata[i] <- image_filenames_grouped[[sample_number]][1]  # Use the first image URL for each sample
  } else {
    customdata[i] <- NA
  }
}

# Add image filenames to the PCA data
PCA_all.df$imagefile <- customdata

# Create a column in PCA data for individual ID (this is so we can join the PCA data to the metadata - some samples had inconsistent sample names)
PCA_all.df <- PCA_all.df |>
  mutate(individualID = gsub(".*_.*_(\\d+).*", "\\1", Sample_name))

# Join both PCA and metadata files together
all <- full_join(Stylo_clusters_metadata, PCA_all.df, by = "individualID")


# -----------------------------------------------------------------------------
# CREATE PLOT -----------------------------------------------------------------
# -----------------------------------------------------------------------------

# dependency code
# Create HTML dependency for the D3 library (version 7.3)
d3 <- htmltools::htmlDependency(
  "d3", "7.3",
  src = c(href = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.3.0/"),  # Specify the source URL for D3 library
  script = "d3.min.js"  # Specify the JavaScript file to load
)

# javascript code 
js <- 'function(el) {
  // Create a tooltip div and assign it the "my-custom-tooltip" class
  var tooltip = d3.select("#" + el.id + " .svg-container")
    .append("div")
    .attr("class", "my-custom-tooltip");

  // Add an event listener for the "plotly_hover" event
  el.on("plotly_hover", function(d) {
    var pt = d.points[0];  // Get the first point of the hover event

    // Set the desired x and y coordinates for the tooltip
    var xPixel = 20; // Change this to 15
    var yPixel = 4;  // Change this to 6

    // Create an image tag with a custom data source and width
    var img = "<img src=\\\"" +  pt.customdata + "\\\" width=400>";
    tooltip.html(img)
      .style("position", "absolute")
      .style("right", xPixel + "px")
      .style("top", yPixel + "px");

    // Apply a transition to the tooltip (so that it fades in)
    tooltip.transition()
      .duration(300)
      .style("opacity", 1);
  });

  // Add an event listener for the "plotly_unhover" event
  el.on("plotly_unhover", function(d) {
    // Apply a transition to hide the tooltip
    tooltip.transition()
      .duration(500)
      .style("opacity", 0);
  });
}'


# create uris (this locates the images so they can be used in the plot)
uris <- purrr::map_chr(
  all$New_sample_name, ~base64enc::dataURI(file = sprintf("Spis_images/resize/%s.JPG", .x))
)


# plot figure
fig <- plot_ly(
  data = all,
  x = ~PC1,
  y = ~PC2,
  color = ~Cluster,
  colors = c("#834177", "#d23359", "#17697c", "#ba4b05", "#00652e", "#b99a2d"),
  customdata = ~uris,
  text = ~paste("Sample Name: ", all$New_sample_name, "<br>",
                "Locality: ", all$locality, "<br>",
                "Site: ", all$EcoLocationID_short, "<br>",
                "Taxa: ", all$Cluster),
  hovertemplate = paste0(
    "<span style='fill:white;font-size:1em;'>%{text}</span><extra></extra>")
) |>
  layout(
    xaxis = list(title = paste0("PC1 (", 
                               signif(pve_all$X100...PCA_all.eig.1.10..sum.PCA_all.eig.1.10..[1],3),"%)"), 
                 zeroline = FALSE, 
                 titlefont = list(size = 15, weight = "bold", color = "black")),
    yaxis = list(title = paste0("PC2 (", 
                                signif(pve_all$X100...PCA_all.eig.1.10..sum.PCA_all.eig.1.10..[2], 3),"%)"), 
                 zeroline = FALSE, 
                 titlefont = list(size = 15, weight = "bold", color = "black")),
    legend = list(x = 15, y = 0.5, color = "black")
  ) |>
  add_markers(x = ~PC1, 
              y = ~PC2, 
              customdata = uris, 
              marker = list(
    size = 13,      # Set the size of the points (adjust as needed)
    opacity = 0.7  # Set the opacity of the points (adjust as needed)
  )) |>
  htmlwidgets::onRender(js)

# add dependencies
fig$dependencies <- c(fig$dependencies, list(d3))

# plot figure
fig

# save fig
htmlwidgets::saveWidget(as_widget(fig), "Spis_PCA_test.html")


