### Below is quarto voodoo to make sure the environment is saved and loaded between pages, because quarto is KM, not MK

#Load common environment if it exists
if (file.exists(here::here("data/env.RData"))) {
  load(here::here("data/env.RData"))
}

#Hook knitr to save common environment on each document render end
local({
  hook_old <- knitr::knit_hooks$get("document") #Grab current hook
  knitr::knit_hooks$set(
    document = function(x) {
      save.image(file = here::here("data/env.RData")) #Save environment
      hook_old(x) #Important: call the original hook, or everything will break
    }
  )
})

### End of quarto voodoo

###Helper lines may be uncommented to install packages manually!

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("org.Hs.eg.db", force = TRUE)
# BiocManager::install("topGO", force = TRUE)
# renv::install("IOR-Bioinformatics/PCSF", force = TRUE)
# renv::install("CogDisResLab/Kinograte", force = TRUE)

###Load all packages statistically for renv
suppressPackageStartupMessages({
  library(here)
  library(babelgene)
  library(tidyverse)
  library(kinograte)
  library(HGNChelper)
  library(creedenzymatic)
  library(DT)
  library(htmltools)
  library(enrichR)
  library(ggnetwork)
  library(igraph)
  library(intergraph)
  library(network)
  library(scatterpie)
  library(fgsea)
  library(PAVER)
  library(ggprism)
  library(svglite)
  library(httr)
  library(jsonlite)
})

knitr::opts_chunk$set(
  comment = "#>",
  cache = FALSE,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

options(readr.show_col_types = FALSE)

set.seed(123)

# Quiet any noisy function
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Make a formatted HTML table
make_table <- function(df, caption=NULL) {
  DT::datatable(
    df,
    rownames = FALSE,
    caption = tags$caption(style = 'text-align: left;', caption),
    options = list(
      scrollX = TRUE,
      scrollY = TRUE,
      paging = TRUE,
      fixedHeader = TRUE,
      pageLength = 10
    )
  ) %>% formatStyle(names(df), "white-space" = "nowrap")
}