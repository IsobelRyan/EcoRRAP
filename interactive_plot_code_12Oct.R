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

Stylo_clusters_metadata <- read_xlsx("Spis_taxa_metadata_13Oct_GRAPH.xlsx")

## add new spreadsheet 

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
directory <- "Spis_images"

# list the files in the directory with the respective pattern
mid <- list.files(directory, pattern = "_MID[0-9]*.JPG", full.names = TRUE)
macro <- list.files(directory, pattern = "_MACRO[0-9]*.JPG", full.names = TRUE)

# remove files with the respective pattern 
file.remove(mid, macro)


# once the colony images are left, now want to remove _COLONY and any digits so that only the sample_name is left 
files <- list.files(directory, full.names = TRUE)


for (file in files) {
  new_name <- gsub("_COLONY[0-9]*", "", basename(file)) # also need to change CBHE to HER
  new_path <- file.path(dirname(file), new_name)
  file.rename(file, new_path)
  cat("Renamed:", file, "to", new_path, "\n")
}

for (file in files) {
  new_name <- gsub("TSAU", "AUK", basename(file)) # also need to change CBHE to HER
  new_path <- file.path(dirname(file), new_name)
  file.rename(file, new_path)
  cat("Renamed:", file, "to", new_path, "\n")
}


all$Sample_name <- gsub("TSDU", "DUN", all$Sample_name)
all$Sample_name <- gsub("ONLI", "LIZ", all$Sample_name)
all$Sample_name <- gsub("TSAU", "AUK", all$Sample_name)
all$Sample_name <- gsub("ONMO", "MOO", all$Sample_name)
all$Sample_name <- gsub("TSMA", "MAS", all$Sample_name)
all$Sample_name <- gsub("CBHE", "HER", all$Sample_name)


# -----------------------------------------------------------------------------
# CREATE COLUMN FOR IMAGE FILE PATHS ------------------------------------------
# -----------------------------------------------------------------------------


# assign Sample_name as row names
PCA_all.df$Sample_name <- rownames(PCA_all.df)

# Rename sample names to the original location names
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

# Function to extract the sample numbers from the filenames
extract_sample_number <- function(filename) {
  number <- sub("^[^0-9]*([0-9]+).*", "\\1", basename(filename))
  return(number)
}

# Extract the sample numbers from the filenames
sample_numbers <- sapply(image_filenames, extract_sample_number)

# Group the image filenames by their sample numbers
image_filenames_grouped <- split(image_filenames, sample_numbers)


customdata <- character(length(PCA_all.df$Sample_name))
for (i in seq_along(PCA_all.df$Sample_name)) {
  sample_number <- extract_sample_number(PCA_all.df$Sample_name[i])
  if (sample_number %in% names(image_filenames_grouped)) {
    customdata[i] <- image_filenames_grouped[[sample_number]][1]  # Use the first image URL for each sample
  } else {
    customdata[i] <- NA
  }
}

PCA_all.df$imagefile <- customdata

PCA_all.df <- PCA_all.df |>
  mutate(individualID = gsub(".*_.*_(\\d+).*", "\\1", Sample_name))


# Join both PCA and metadata files together
all <- full_join(Stylo_clusters_metadata, PCA_all.df, by = "individualID")

# -----------------------------------------------------------------------------
# FOR RESIZING IMAGES ---------------------------------------------------------
# -----------------------------------------------------------------------------

for (i in 284:nrow(PCA_all.df)) {
  # Get the file path from the DataFrame
  file_path <- PCA_all.df$imagefile[i]
  
  # Load the image
  img <- image_read(file_path)
  
  # Resize the image to 50% with Lanczos filter (adjust as needed)
  resized_img <- image_resize(img, geometry = "20%")
  
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Define the path to save the resized image
  save_path <- paste0("Spis_images/resize/", file_name, ".JPG")
  
  # Save the resized image
  image_write(resized_img, path = save_path)
  
  # Print a message to indicate progress
  cat("Processed:", file_path, "\n")
}


# -----------------------------------------------------------------------------
# CREATE PLOT -----------------------------------------------------------------
# -----------------------------------------------------------------------------

# dependency 
d3 <- htmltools::htmlDependency(
  "d3", "7.3",
  src = c(href = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.3.0/"),
  script = "d3.min.js"
)


# javascript code 
js <- 'function(el) {
  var tooltip = d3.select("#" + el.id + " .svg-container")
  .append("div")
  .attr("class", "my-custom-tooltip");
  
  el.on("plotly_hover", function(d) {
    var pt = d.points[0];
    
    // Set the desired x and y coordinates
    var xPixel = 20; // Change this to 15
    var yPixel = 4;  // Change this to 6
    
    var img = "<img src=\\\"" +  pt.customdata + "\\\" width=400>";
    tooltip.html(img)
    .style("position", "absolute")
    .style("right", xPixel + "px")
    .style("top", yPixel + "px");
    tooltip.transition()
    .duration(300)
    .style("opacity", 1);
  });
  
  el.on("plotly_unhover", function(d) {
    tooltip.transition()
    .duration(500)
    .style("opacity", 0);
  });
}'





# create uris
uris <- purrr::map_chr(
  all$New_sample_name, ~base64enc::dataURI(file = sprintf("https://raw.githubusercontent.com/IsobelRyan/EcoRRAP/main/Spis_images/%s.JPG", .x))
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

fig$dependencies <- c(fig$dependencies, list(d3))

fig

# save fig
htmlwidgets::saveWidget(as_widget(fig), "Spis_PCA_Interactive_Plot.html")

write_html

all <- all[1:379, ]
