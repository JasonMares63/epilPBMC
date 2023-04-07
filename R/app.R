
#############
#setwd("Epilepsy_PBMC")

runShiny <- function(){


require(tidyverse)
require(stringr)
require(RColorBrewer)
require(shiny)
require(patchwork)
require(Seurat)
require(DT)
require(ggrepel)
require(dplyr)
require(ggpubr)
require(ggprism)
require(shinythemes)

source("R/global.R")
source("R/ui.R")
source("R/server.R")

shinyApp(ui, server)

}
