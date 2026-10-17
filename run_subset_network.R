source("subset_network.R")

#Gene sets specified from grant fig
gene_sets <- list(
  "Actin Turnover" = c(paste0("PAK", 1:4), paste0("LIMK", 1:2), paste0("CFL", 1:2)),
  "Actin Branching / Spine Structure" = c(paste0("CYFIP", 1:2), paste0("ABI", 1:3), paste0("WASF", 1:3), "ACTR2", "ACTR3"),
  "Immune / ROS Signaling" = c("NCF1", "NCF2", "NCF4", "CYBB")
)

scores <- read.csv("results/combined_node_scores.csv", stringsAsFactors = FALSE)
is_kinase <- vapply(strsplit(ifelse(is.na(scores$type), "", scores$type), "|", fixed = TRUE),
                    function(v) any(tolower(trimws(v)) == "kinase"), logical(1))
kinases <- unique(scores$name[is_kinase & !is.na(scores$name) & nzchar(trimws(scores$name))])
groups <- c(gene_sets, list(Kinase = kinases))

x <- plot_gene_set(
  gene_sets = gene_sets,
  remove_isolates = FALSE,
  expand_one_hop = TRUE,
  groups = groups,
  layout_by_group = TRUE,
  group_layout = c(
    "Actin Turnover" = "spatial",
    "Actin Branching / Spine Structure" = "spatial",
    "Immune / ROS Signaling" = "spatial"
    #omitted kinase group so its mode defaults to "connectivity"
  ),
  color_node_groups = FALSE,
  node_groups = groups,
  color_by_annotation = TRUE,
  group_colors = c(
    "Actin Turnover" = "#0072B2",
    "Actin Branching / Spine Structure" = "#009E73",
    "Immune / ROS Signaling" = "#CC79A7",
    "Kinase" = "#D55E00"
  ),
  title = "",
  out = "results/hypothesis_driven2"
)
print(x$plot)

go <- "GO:0016601"
db <- msigdbr::msigdbr(collection = "C5")
RAC1_genes <- unique(db$gene_symbol[db$gs_exact_source == go])
go <- "GO:0007015"
actin_genes <- unique(db$gene_symbol[db$gs_exact_source == go])
go <- "GO:0014069"
synaptic_genes <- unique(db$gene_symbol[db$gs_exact_source == go])
go <- "GO:0034614"
ROS_genes <- unique(db$gene_symbol[db$gs_exact_source == go])
go <- "GO:0048167"
synaptic_genes <- unique(db$gene_symbol[db$gs_exact_source == go])

gene_sets <- list(
  #"Rac Protein Signaling Transduction" = RAC1_genes,
  #"Actin Filament Organization" = actin_genes,
  "Regulation of Synaptic Plasticity" = synaptic_genes
  #"Cellular Response to ROS" = ROS_genes
)

groups <- c(gene_sets)#, list(Kinase = kinases))

x <- plot_gene_set(
  gene_sets = gene_sets,
  remove_isolates = TRUE,
  expand_one_hop = FALSE,
  groups = groups,
  layout_by_group = TRUE,
  group_layout = c(
    #"Rac Protein Signaling Transduction" = 'spatial',
    #"Actin Filament Organization" = 'spatial',
    "Regulation of Synaptic Plasticity" = 'spatial'
    #"Cellular Response to ROS" = 'spatial'
    # Kinase is omitted, so its mode defaults to "connectivity".
  ),
  color_node_groups = FALSE,
  node_groups = groups,
  color_by_annotation = TRUE,
  group_colors = c(
    #"Rac Protein Signaling Transduction" = "#0072B2",
    #"Actin Filament Organization" = "#009E73",
    "Regulation of Synaptic Plasticity" = "#CC79A7"
    #"Cellular Response to ROS" = "#A64452"
    #"Kinase" = "#D55E00"
  ),
  title = "",
  out = "results/gene_sets"
)
print(x$plot)

#Initial list w one-hop neighbors added
gene_sets <- list(
  "Actin Turnover" = c(paste0("PAK", 1:4), paste0("LIMK", 1:2), paste0("CFL", 1:2), paste0("RAC", 1:3), 
                       paste0("BMPR2"), paste0("ARHGEF7"), paste0("CDC42") ,paste0("NCK1"), paste0("TAGLN2"),paste0("ROCK1")),
  "Actin Branching / Spine Structure" = c(paste0("CYFIP", 1:2), paste0("ABI", 1:3), paste0("WASF", 1:3), "ACTR2", "ACTR3",
                                          "ABL1", "PCDH10", "HCSL1", "CORO1C"),
  "Immune / ROS Signaling" = c("NCF1", "NCF2", "NCF4", "CYBB"),
  "Regulation of Synaptic Plasticity" = synaptic_genes
)

groups <- c(gene_sets)#, list(Kinase = kinases))

x <- plot_gene_set(
  gene_sets = gene_sets,
  remove_isolates = TRUE,
  expand_one_hop = FALSE,
  groups = groups,
  layout_by_group = TRUE,
  group_layout = c(
    "Actin Turnover" = "spatial",
    "Actin Branching / Spine Structure" = "spatial",
    "Immune / ROS Signaling" = "spatial",
    "Regulation of Synaptic Plasticity" = 'spatial'
  ),
  color_node_groups = FALSE,
  node_groups = groups,
  color_by_annotation = TRUE,
  group_colors = c(
    "Actin Turnover" = "#0072B2",
    "Actin Branching / Spine Structure" = "#009E73",
    "Immune / ROS Signaling" = "#CC79A7",
    "Regulation of Synaptic Plasticity" = "#A64452"
    #"Kinase" = "#D55E00"
  ),
  title = "",
  out = "results/john_request"
)
print(x$plot)
