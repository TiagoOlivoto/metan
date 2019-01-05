## ---- eval = FALSE, collapse=TRUE, comment = "#"-------------------------
#  # download the package from Github
#  devtools::install_github("TiagoOlivoto/METAAB")
#  

## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
require(METAAB)
require(ggplot2)
require(cowplot) # used in this material to arrange the graphics
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/2/files/1561de7f-b8fd-4093-b4d1-bfcef299dd22/WAASBdata.csv?dl=1")


