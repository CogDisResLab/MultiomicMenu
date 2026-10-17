# Standalone: needs igraph, ggplot2 and ggrepel. No Quarto objects required.
# source() defines functions only; see run_subset_network.R for an example.
# groups controls layout; node_groups/group_colors control color.
# With layout_by_group = TRUE, group_layout defaults to "connectivity".
# Use group_layout = c(MyGroup = "spatial") to cluster selected groups.
REMOVE_ISOLATES <- FALSE
EXPAND_ONE_HOP <- FALSE
COLOR_NODE_GROUPS <- FALSE
COLOR_BY_ANNOTATION <- FALSE
LAYOUT_BY_GROUP <- FALSE
SHOW_LEGEND <- TRUE

# source("subset_network.R"); x <- subset_network(); print(x$plot)
# Normalize a gene vector or named list without relying on any session objects.
.network_sets <- function(x, allow_reserved = FALSE) {
  if (is.character(x)) x <- list(Selected = x)
  if (!is.list(x) || !length(x) || is.null(names(x)))
    stop("gene_sets must be a gene vector or a named list of gene vectors")
  names(x) <- trimws(names(x))
  if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
      (!allow_reserved && any(names(x) %in% c("Shared", "Neighbor"))))
    stop("Gene-set names must be unique, nonblank, and not Shared or Neighbor")
  x <- lapply(x, function(v) {
    if (!is.character(v)) stop("Each gene set must be a character vector of gene names")
    v <- toupper(trimws(v)); unique(v[!is.na(v) & nzchar(v)])
  })
  if (any(!lengths(x))) stop("Gene sets must not be empty")
  x
}

# Join annotation tokens without duplicating vertices or changing structural node type.
.network_annotations <- function(node_names, file) {
  if (!is.character(file) || length(file) != 1L || is.na(file))
    stop("annotations_file must be a single file path")
  if (!file.exists(file)) stop("Annotation file not found: ", file)
  scores <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("name", "type") %in% names(scores)))
    stop("Annotation file must contain name and type columns")
  key <- function(z) toupper(trimws(as.character(z)))
  ids <- key(scores$name)
  if (anyNA(ids) || any(!nzchar(ids))) stop("Annotation names must be nonblank")
  rows <- split(seq_len(nrow(scores)), ids)
  known <- c("Protein", "Kinase", "TF", "RNA", "Peptide", "Phosphosite", "Hidden")
  sets <- lapply(key(node_names), function(id) {
    types <- as.character(scores$type[rows[[id]]])
    tokens <- trimws(unlist(strsplit(types[!is.na(types)], "|", fixed = TRUE)))
    tokens <- tokens[nzchar(tokens) & toupper(tokens) != "NA"]
    matched <- match(toupper(tokens), toupper(known))
    tokens[!is.na(matched)] <- known[matched[!is.na(matched)]]
    tokens <- unique(tokens)
    if (length(tokens)) tokens else "Hidden"
  })
  palette <- c(Protein = "#E69F00", Kinase = "#0072B2", TF = "#CC79A7", RNA = "#009E73",
               Peptide = "#56B4E9", Phosphosite = "#D55E00", Hidden = "#BDBDBD")
  extra <- sort(setdiff(unique(unlist(sets)), names(palette)))
  if (length(extra)) palette <- c(palette, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  list(sets = sets, colors = palette)
}

# Unspecified groups stay connectivity-based; named entries override that default.
.network_layout_modes <- function(group_names, group_layout) {
  modes <- setNames(rep("connectivity", length(group_names)), group_names)
  if (is.null(group_layout)) return(modes)
  if (!is.character(group_layout) || !length(group_layout) || anyNA(group_layout) ||
      any(!group_layout %in% c("connectivity", "spatial")))
    stop('group_layout must contain "connectivity" or "spatial"')
  if (is.null(names(group_layout))) {
    if (length(group_layout) != 1L) stop("Unnamed group_layout must be a single mode")
    modes[] <- group_layout
  } else {
    ids <- names(group_layout)
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("group_layout names must be unique and nonblank")
    if (any(!ids %in% group_names))
      stop("Unknown group_layout groups: ", paste(setdiff(ids, group_names), collapse = ", "))
    modes[ids] <- group_layout
  }
  modes
}

# Relax free nodes around spatial anchors using the original, unweighted topology.
.network_connectivity_layout <- function(g, xy, anchored, base_xy) {
  if (!any(anchored) || all(anchored) || !igraph::ecount(g)) return(xy)
  # Disconnected components also remain free to avoid collisions with moved anchors.
  fixed <- anchored
  ends <- igraph::ends(g, igraph::E(g), names = FALSE)
  lengths <- sqrt(rowSums((base_xy[ends[, 1], , drop = FALSE] -
                           base_xy[ends[, 2], , drop = FALSE])^2))
  unit <- stats::median(lengths[is.finite(lengths) & lengths > 0])
  if (!is.finite(unit) || unit <= 0) unit <- 1
  start <- xy / unit
  lower_x <- lower_y <- rep(-Inf, nrow(xy))
  upper_x <- upper_y <- rep(Inf, nrow(xy))
  lower_x[fixed] <- upper_x[fixed] <- start[fixed, 1]
  lower_y[fixed] <- upper_y[fixed] <- start[fixed, 2]
  result <- igraph::layout_with_fr(g, coords = start, weights = NA,
    niter = 1500, grid = "nogrid", minx = lower_x, maxx = upper_x,
    miny = lower_y, maxy = upper_y) * unit
  result[fixed, ] <- xy[fixed, , drop = FALSE]
  if (any(!is.finite(result))) stop("Grouped connectivity layout produced non-finite coordinates")
  result
}

# Place explicit spatial memberships independently of the displayed colors.
.network_group_layout <- function(xy, membership) {
  membership <- membership[, colSums(membership) > 0, drop = FALSE]
  groups <- colnames(membership)
  centers <- data.frame(group = character(), x = numeric(), y = numeric())
  if (!length(groups)) return(list(xy = xy, centers = centers))
  angles <- pi / 2 + 2 * pi * (seq_along(groups) - 1) / length(groups)
  centers <- data.frame(group = groups, x = 0.78 * cos(angles), y = 0.78 * sin(angles))
  if (length(groups) == 2L) {
    centers$x <- c(-0.78, 0.78); centers$y <- 0
  }
  if (length(groups) == 1L) centers[c("x", "y")] <- 0
  anchored <- rowSums(membership) > 0
  mask <- membership[anchored, , drop = FALSE]
  combinations <- apply(mask, 1, paste0, collapse = "")
  # Give larger groups more room while keeping the existing node radius unchanged.
  largest_spread <- max(0.22, 0.055 * sqrt(max(table(combinations))))
  center_scale <- max(1, 3 * largest_spread / 0.78)
  centers[c("x", "y")] <- centers[c("x", "y")] * center_scale
  targets <- (mask / rowSums(mask)) %*% as.matrix(centers[c("x", "y")])
  result <- targets
  source_xy <- xy[anchored, , drop = FALSE]
  for (ids in split(seq_len(nrow(source_xy)), combinations)) {
    if (length(ids) < 2L) next
    local <- sweep(source_xy[ids, , drop = FALSE], 2, colMeans(source_xy[ids, , drop = FALSE]))
    angular_order <- order(atan2(local[, 2], local[, 1]), ids)
    count <- length(ids)
    theta <- if (count <= 8L) 2 * pi * (seq_len(count) - 1) / count else
      (seq_len(count) - 1) * pi * (3 - sqrt(5))
    radius <- if (count <= 8L) rep(1, count) else sqrt((seq_len(count) - 0.5) / count)
    local[angular_order, ] <- cbind(radius * cos(theta), radius * sin(theta))
    local <- sweep(local, 2, colMeans(local))
    local <- local / max(sqrt(rowSums(local^2)))
    spread <- max(0.45, 0.090 * sqrt(count))
    result[ids, ] <- targets[ids, , drop = FALSE] + spread * local
  }
  xy[anchored, ] <- result
  list(xy = xy, centers = centers)
}

# Small geometry helper: equal sectors encode membership, not quantities.
.network_pies <- function(n, plot_sets) {
  slices <- circles <- vector("list", nrow(n))
  for (i in seq_len(nrow(n))) {
    count <- length(plot_sets[[i]])
    slices[[i]] <- do.call(rbind, lapply(seq_len(count), function(j) {
      theta <- seq(pi / 2 + 2 * pi * (j - 1) / count,
                   pi / 2 + 2 * pi * j / count, length.out = 65)
      px <- n$x[i] + n$plot_radius[i] * cos(theta)
      py <- n$y[i] + n$plot_radius[i] * sin(theta)
      if (count > 1L) {
        px <- c(n$x[i], px, n$x[i]); py <- c(n$y[i], py, n$y[i])
      }
      data.frame(x = px, y = py,
                 node = n$name[i], piece = paste(i, j, sep = ":"),
                 group = plot_sets[[i]][j])
    }))
    theta <- seq(0, 2 * pi, length.out = 97)
    circles[[i]] <- data.frame(x = n$x[i] + n$plot_radius[i] * cos(theta),
                              y = n$y[i] + n$plot_radius[i] * sin(theta), node = n$name[i])
  }
  list(slices = do.call(rbind, slices), circles = do.call(rbind, circles))
}

# Plotting stays independent of selection, expansion and annotation logic.
.network_plot <- function(n, e, plot_sets, colors, title, subtitle, caption,
                          seed, show_legend, font_family) {
  pies <- .network_pies(n, plot_sets)
  used <- names(colors)[names(colors) %in% unique(unlist(plot_sets))]
  # Self-loops have no distinct endpoint and cannot be drawn by geom_curve.
  drawn_edges <- e[e$from != e$to, , drop = FALSE]
  limits <- function(z) {
    half <- max(diff(range(z)) / 2 + 0.18, 0.45)
    mean(range(z)) + c(-half, half)
  }
  ggplot2::ggplot() +
    ggplot2::geom_curve(data = drawn_edges, ggplot2::aes(x, y, xend = xend, yend = yend),
      curvature = 0.05, colour = "grey50", alpha = 0.5, linewidth = 0.35) +
    ggplot2::geom_polygon(data = pies$slices, ggplot2::aes(x, y, group = piece, fill = group),
      colour = "white", linewidth = 0.12, show.legend = FALSE) +
    ggplot2::geom_path(data = pies$circles, ggplot2::aes(x, y, group = node),
      colour = "#444444", linewidth = 0.25) +
    # Invisible points provide clean circular legend keys for the pie colors.
    ggplot2::geom_point(data = n, ggplot2::aes(x, y, fill = group), shape = 21, size = 0, alpha = 0,
      colour = "#444444", show.legend = TRUE) +
    ggrepel::geom_text_repel(data = n, ggplot2::aes(x, y, label = name),
      family = font_family, fontface = "bold.italic", size = 2.8, seed = seed,
      box.padding = 0.3, point.padding = 0.55, point.size = 2.5, force = 1.8,
      max.overlaps = Inf, max.iter = 20000, max.time = 2,
      min.segment.length = 0.1, segment.colour = "grey65", segment.size = 0.25) +
    ggplot2::scale_fill_manual(values = colors, limits = used, name = NULL,
      labels = function(z) vapply(z, function(v) paste(strwrap(v, 24), collapse = "\n"), "")) +
    ggplot2::coord_equal(xlim = limits(n$x), ylim = limits(n$y), clip = "off", expand = FALSE) +
    ggplot2::theme_void(base_size = 10, base_family = font_family) +
    ggplot2::labs(title = paste(strwrap(title, 60), collapse = "\n"), subtitle = NULL,
      caption = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9, margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(colour = "grey40", size = 8, hjust = 0, margin = ggplot2::margin(t = 8)),
      plot.margin = ggplot2::margin(16, 18, 12, 18),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = ggplot2::element_text(size = 9), legend.key.width = grid::unit(0.55, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = ceiling(length(used) / 4), byrow = TRUE,
      override.aes = list(shape = 21, size = 3.6, alpha = 1)))
}

subset_network <- function(
    results = "results", enrich_file = file.path(results, "enrichr_GO.csv"),
    mode = c("calcium_signaling", "calcium", "neuronal", "custom"),
    pathways = NULL, gmt = NULL, go_only = TRUE, seed = 123,
    out = file.path(results, "subset_network"), title = NULL,
    remove_isolates = REMOVE_ISOLATES, expand_one_hop = EXPAND_ONE_HOP,
    gene_sets = NULL, color_node_groups = COLOR_NODE_GROUPS,
    node_groups = NULL, group_colors = NULL, show_legend = SHOW_LEGEND,
    font_family = "sans", width = 8, height = 6,
    color_by_annotation = COLOR_BY_ANNOTATION,
    annotations_file = file.path(results, "combined_node_scores.csv"),
    layout_by_group = LAYOUT_BY_GROUP, groups = NULL,
    group_layout = "connectivity") {
  mode <- match.arg(mode)
  if (!is.logical(remove_isolates) || length(remove_isolates) != 1L || is.na(remove_isolates))
    stop("remove_isolates must be TRUE or FALSE")
  if (!is.logical(expand_one_hop) || length(expand_one_hop) != 1L || is.na(expand_one_hop))
    stop("expand_one_hop must be TRUE or FALSE")
  if (!is.logical(color_node_groups) || length(color_node_groups) != 1L || is.na(color_node_groups))
    stop("color_node_groups must be TRUE or FALSE")
  if (!is.logical(show_legend) || length(show_legend) != 1L || is.na(show_legend))
    stop("show_legend must be TRUE or FALSE")
  if (!is.logical(go_only) || length(go_only) != 1L || is.na(go_only))
    stop("go_only must be TRUE or FALSE")
  if (!is.character(font_family) || length(font_family) != 1L || is.na(font_family))
    stop("font_family must be one font-family name")
  for (size in list(width, height))
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0)
      stop("width and height must be positive finite numbers")
  if (!is.logical(color_by_annotation) || length(color_by_annotation) != 1L || is.na(color_by_annotation))
    stop("color_by_annotation must be TRUE or FALSE")
  if (!is.logical(layout_by_group) || length(layout_by_group) != 1L || is.na(layout_by_group))
    stop("layout_by_group must be TRUE or FALSE")
  if (color_node_groups && !color_by_annotation && is.null(node_groups))
    stop("Provide node_groups when color_node_groups is TRUE")
  for (p in c("igraph", "ggplot2", "ggrepel"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
  read <- function(f) {
    if (!file.exists(f)) stop("File not found: ", f)
    utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  }
  col <- function(d, choices) {
    clean <- function(z) gsub("[^a-z0-9]", "", tolower(z))
    i <- match(clean(choices), clean(names(d))); i <- i[!is.na(i)]
    if (!length(i)) stop("Missing column: ", paste(choices, collapse = " / "))
    d[[i[1]]]
  }
  key <- function(z) toupper(trimws(as.character(z)))
  split_genes <- function(z) {
    z <- key(unlist(strsplit(ifelse(is.na(z), "", z), "[;,|[:space:]]+")))
    unique(z[nzchar(z) & z != "NA"])
  }
  n <- read(file.path(results, "network_nodes.csv"))
  e <- read(file.path(results, "network_edges.csv"))
  if (!nrow(n)) stop("Node table is empty")
  if (!"name" %in% names(n) || !all(c("from", "to") %in% names(e)))
    stop("Expected nodes: name; edges: from, to")
  n <- n[c("name", setdiff(names(n), "name"))]
  e <- e[c("from", "to", setdiff(names(e), c("from", "to")))]
  n$name <- trimws(as.character(n$name))
  e$from <- trimws(as.character(e$from)); e$to <- trimws(as.character(e$to))
  if (anyNA(n$name) || any(!nzchar(n$name)) || anyDuplicated(key(n$name)))
    stop("Node names must be nonblank and unique (case-insensitive)")
  if (anyNA(e$from) || anyNA(e$to) ||
      !all(c(e$from, e$to) %in% n$name)) stop("Invalid/missing edge endpoints")
  if (!"type" %in% names(n)) n$type <- "Node"
  n$type[is.na(n$type) | !nzchar(n$type)] <- "Node"
  edge_columns <- names(e)

  direct_selection <- !is.null(gene_sets)
  if (direct_selection) {
    sets <- .network_sets(gene_sets)
    selected <- data.frame()
  } else {
    a <- read(enrich_file)
    a$.term <- as.character(col(a, c("Term", "Description", "Pathway")))
    raw_fdr <- col(a, c("Adjusted.P.value", "Adjusted P-value", "adjusted_pvalue", "FDR", "padj", "qvalue"))
    a$.fdr <- suppressWarnings(as.numeric(raw_fdr))
    if (any(!is.na(raw_fdr) & nzchar(trimws(as.character(raw_fdr))) & is.na(a$.fdr)))
      stop("Non-numeric adjusted P values; use numeric Enrichr exports")
    if (any(!is.na(a$.fdr) & (!is.finite(a$.fdr) | a$.fdr < 0 | a$.fdr > 1)))
      stop("Adjusted P values must be numeric in [0,1]")
    if (go_only) {
      is_go <- grepl("GO:[0-9]{7}", a$.term)
      if ("Library" %in% names(a)) is_go <- is_go | grepl("^GO_", a$Library)
    } else is_go <- rep(TRUE, nrow(a))
    a <- a[!is.na(a$.term) & !is.na(a$.fdr) & a$.fdr < 0.05 &
             !is.na(is_go) & is_go, , drop = FALSE]
    if (!nrow(a)) stop("No qualifying GO terms with adjusted P < 0.05")
    term <- tolower(gsub("[-_]", " ", a$.term))
    ca <- grepl("calcium|calmodulin|calcineurin|calbindin|(^|[^a-z])ca2?\\+", term)
    signal <- grepl("signal|signalling|messenger", term)
    neuro <- grepl("neuron|neural|synap|neurotransmi|neuropeptide|neurotrophin", term) &
      grepl("signal|signalling|transmi|synaptic (vesicle|plasticity)|postsynaptic|presynaptic", term)
    if (mode == "custom") {
      if (is.null(pathways) || !length(pathways) || is.null(names(pathways)) ||
          anyNA(names(pathways)) || any(!nzchar(names(pathways))) ||
          any(names(pathways) %in% c("Shared", "Neighbor")) || anyDuplicated(names(pathways)))
        stop("custom pathways must be a named list of GO IDs or exact term names")
      term_groups <- lapply(pathways, function(ids) {
        goid <- ifelse(grepl("GO:[0-9]{7}", a$.term),
                       sub(".*(GO:[0-9]{7}).*", "\\1", a$.term), "")
        which(key(a$.term) %in% key(ids) | key(goid) %in% key(ids))
      })
      if (any(!lengths(term_groups))) stop("No significant matching terms for: ",
                                     paste(names(term_groups)[!lengths(term_groups)], collapse = ", "))
    } else {
      hit <- switch(mode, calcium = ca, calcium_signaling = ca & signal, neuronal = neuro)
      term_groups <- setNames(list(which(hit)), gsub("_", " ", mode))
      if (!any(hit)) stop("No significant terms match mode: ", mode)
    }
    # Enrichr Genes are QUERY OVERLAPS. For a different query, pass the matching
    # Enrichr GMT so every whole-network member of a significant term is eligible.
    if (is.null(gmt)) {
      a$.genes <- as.character(col(a, c("Genes", "Overlapping Genes", "overlapping_genes")))
      members <- lapply(a$.genes, split_genes)
    } else {
      if (!file.exists(gmt)) stop("GMT not found: ", gmt)
      z <- strsplit(readLines(gmt, warn = FALSE), "\t", fixed = TRUE)
      z <- z[lengths(z) >= 3]
      terms <- vapply(z, `[`, "", 1)
      term_key <- function(t) ifelse(grepl("GO:[0-9]{7}", t),
        sub(".*(GO:[0-9]{7}).*", "\\1", t), key(t))
      ix <- match(term_key(a$.term), term_key(terms))
      wanted <- unique(unlist(term_groups))
      if (anyNA(ix[wanted])) stop("Selected term missing from GMT; check library/version")
      members <- lapply(ix, function(i) if (is.na(i)) character() else split_genes(z[[i]][-c(1,2)]))
    }
    sets <- lapply(term_groups, function(i) unique(unlist(members[i])))
    selected_indices <- sort(unique(unlist(term_groups)))
    selected <- a[selected_indices, , drop = FALSE]
    selected$matched_groups <- vapply(selected_indices, function(i) {
      paste(names(term_groups)[vapply(term_groups, function(v) i %in% v,
                                logical(1))], collapse = " | ")
    }, "")
  }
  membership <- matrix(unlist(lapply(sets, function(s) key(n$name) %in% s)),
                       nrow = nrow(n), dimnames = list(NULL, names(sets)))
  missing_genes <- lapply(sets, function(v) setdiff(v, key(n$name)))
  if (direct_selection && any(lengths(missing_genes))) {
    missing <- unique(unlist(missing_genes))
    warning(length(missing), " selected genes absent from network: ",
      paste(utils::head(missing, 10), collapse = ", "),
      if (length(missing) > 10) " ... (all names in missing_genes)", call. = FALSE)
  }
  is_seed <- rowSums(membership) > 0
  if (!any(is_seed)) stop("Selected terms have no genes in the network; check identifiers/species")
  keep <- is_seed
  if (expand_one_hop) {
    # Read neighbors from the FULL network, only once, before dropping isolates.
    seeds <- n$name[is_seed]
    neighbors <- unique(c(e$to[e$from %in% seeds], e$from[e$to %in% seeds]))
    keep <- keep | n$name %in% neighbors
  }
  n$is_seed <- is_seed
  n <- n[keep, , drop = FALSE]; membership <- membership[keep, , drop = FALSE]
  n$pathways <- apply(membership, 1, function(v) paste(names(sets)[v], collapse = " | "))
  n$group <- ifelse(!n$is_seed, "Neighbor",
                    ifelse(rowSums(membership) > 1, "Shared", n$pathways))
  e <- e[e$from %in% n$name & e$to %in% n$name, , drop = FALSE]
  # Undirected context edges, no arrows or causal claims. Preserve original metadata.
  g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = n)
  n$degree <- as.numeric(igraph::degree(igraph::simplify(g), loops = FALSE))
  if (remove_isolates) {
    keep_connected <- n$degree > 0
    if (!any(keep_connected))
      stop("No nodes remain after removing isolates; disable remove_isolates or enable one-hop expansion")
    n <- n[keep_connected, , drop = FALSE]
    membership <- membership[keep_connected, , drop = FALSE]
    e <- e[e$from %in% n$name & e$to %in% n$name, , drop = FALSE]
    g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = n)
  }
  set.seed(seed)
  xy <- if (nrow(n) == 1) matrix(c(0, 0), 1) else if (!nrow(e))
    igraph::layout_in_circle(g) else igraph::layout_with_fr(g, weights = NA, niter = 1500)
  # Rotate the longest layout axis horizontally, then normalize without distortion.
  if (nrow(n) > 1) {
    xy <- sweep(xy, 2, colMeans(xy))
    xy <- xy %*% eigen(crossprod(xy), symmetric = TRUE)$vectors
    xy <- sweep(xy, 2, apply(xy, 2, function(v) mean(range(v))))
    span <- max(apply(xy, 2, function(v) diff(range(v))))
    if (span > 0) xy <- 2 * xy / span
  }
  n$x <- xy[, 1]; n$y <- xy[, 2]
  n$plot_radius <- 0.022 + 0.013 * sqrt(n$degree / max(c(1, n$degree)))
  palette <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#333333")
  if (length(sets) > length(palette)) palette <- grDevices::hcl.colors(length(sets), "Dark 3")
  colors <- setNames(palette[seq_along(sets)], names(sets))
  colors <- c(colors, Shared = "#524A67", Neighbor = "#D9D9D9")
  n$override_group <- ""
  if (color_node_groups && !color_by_annotation) {
    overlays <- .network_sets(node_groups)
    extra <- setdiff(names(overlays), names(colors))
    if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), "Set 2"), extra))
    # First group wins if a node appears in more than one override group.
    for (group in rev(names(overlays))) {
      hit <- key(n$name) %in% overlays[[group]]
      n$group[hit] <- group; n$override_group[hit] <- group
    }
  }
  plot_sets <- lapply(seq_len(nrow(n)), function(i) {
    if (nzchar(n$override_group[i])) n$override_group[i] else if (!n$is_seed[i]) "Neighbor" else
      names(sets)[membership[i, ]]
  })
  if (color_by_annotation) {
    annotation <- .network_annotations(n$name, annotations_file)
    plot_sets <- annotation$sets
    colors <- annotation$colors
    n$annotation_type <- vapply(plot_sets, paste, "", collapse = " | ")
    n$group <- vapply(plot_sets, function(v) if (length(v) == 1L) v else "Shared", "")
  }
  if (!is.null(group_colors)) {
    if (!is.character(group_colors) || is.null(names(group_colors)) || anyNA(names(group_colors)) ||
        any(!nzchar(names(group_colors))) || anyDuplicated(names(group_colors)) || anyNA(group_colors))
      stop("group_colors must be a named character vector of valid colors")
    tryCatch(grDevices::col2rgb(group_colors), error = function(e) stop("Invalid group_colors color"))
    colors[names(group_colors)] <- group_colors
  }
  layout_centers <- data.frame(group = character(), x = numeric(), y = numeric())
  layout_membership <- matrix(FALSE, nrow(n), 0, dimnames = list(n$name, NULL))
  layout_group_modes <- setNames(character(), character())
  n$layout_groups <- ""
  n$layout_mode <- "connectivity"
  if (layout_by_group) {
    spatial_groups <- if (is.null(groups)) sets else .network_sets(groups, allow_reserved = TRUE)
    layout_membership <- vapply(spatial_groups, function(v) key(n$name) %in% v,
                                 logical(nrow(n)))
    dim(layout_membership) <- c(nrow(n), length(spatial_groups))
    dimnames(layout_membership) <- list(n$name, names(spatial_groups))
    n$layout_groups <- apply(layout_membership, 1, function(v)
      paste(names(spatial_groups)[v], collapse = " | "))
    layout_group_modes <- .network_layout_modes(names(spatial_groups), group_layout)
    spatial_membership <- layout_membership[, layout_group_modes == "spatial", drop = FALSE]
    anchored <- rowSums(spatial_membership) > 0
    n$layout_mode[anchored] <- "spatial"
    base_xy <- as.matrix(n[c("x", "y")])
    grouped <- .network_group_layout(base_xy, spatial_membership)
    set.seed(seed)
    grouped$xy <- .network_connectivity_layout(g, grouped$xy, anchored, base_xy)
    n$x <- grouped$xy[, 1]; n$y <- grouped$xy[, 2]
    layout_centers <- grouped$centers
  }
  e$x <- n$x[match(e$from, n$name)]; e$y <- n$y[match(e$from, n$name)]
  e$xend <- n$x[match(e$to, n$name)]; e$yend <- n$y[match(e$to, n$name)]
  n$plot_groups <- vapply(plot_sets, paste, "", collapse = " | ")
  g <- igraph::graph_from_data_frame(e[edge_columns], directed = FALSE, vertices = n)
  plot_title <- if (!is.null(title)) title else if (direct_selection) "Selected gene sets" else gsub("_", " ", mode)
  subtitle <- sprintf("%s nodes | %s original edges | %s %s", nrow(n), nrow(e),
    if (direct_selection) length(sets) else nrow(selected),
    if (direct_selection) "gene sets" else "significant terms")
  caption <- paste(if (direct_selection) "Selected gene sets; no enrichment test;" else
    "Enrichr adjusted P < 0.05;", if (expand_one_hop) "one-hop neighbors included;" else "induced subgraph;",
    if (remove_isolates) "isolates removed." else "isolates retained.",
    if (color_by_annotation) "Pie slices show combined-node annotations;" else "Pie slices show membership;",
    "node size indicates subset degree.",
    if (any(n$layout_mode == "spatial")) "Specified spatial groups; other nodes positioned by connectivity." else "")
  p <- .network_plot(n, e, plot_sets, colors, plot_title, subtitle, caption,
                     seed, show_legend, font_family)
  if (!is.null(out)) {
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(paste0(out, ".png"), p, width = width, height = height, dpi = 400, bg = "white")
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(paste0(out, ".pdf"), p, device = pdf_device,
                    width = width, height = height, bg = "white")
    utils::write.csv(n, paste0(out, "_nodes.csv"), row.names = FALSE)
    utils::write.csv(e[edge_columns],
                     paste0(out, "_edges.csv"), row.names = FALSE)
    if (!direct_selection)
      utils::write.csv(selected[!grepl("^\\.", names(selected))], paste0(out, "_terms.csv"), row.names = FALSE)
  }
  invisible(list(plot = p, graph = g, nodes = n, edges = e, terms = if (direct_selection) NULL else selected,
                 membership = membership, missing_genes = missing_genes, colors = colors,
                 plot_sets = plot_sets, edge_columns = edge_columns, layout_centers = layout_centers,
                 layout_membership = layout_membership, layout_group_modes = layout_group_modes))
}

plot_gene_set <- function(
    genes = NULL, results = "results", out = file.path(results, "selected_gene_set"),
    title = "Selected gene set", remove_isolates = REMOVE_ISOLATES,
    expand_one_hop = EXPAND_ONE_HOP, gene_sets = NULL,
    color_node_groups = COLOR_NODE_GROUPS, node_groups = NULL, group_colors = NULL,
    show_legend = SHOW_LEGEND, font_family = "sans", width = 8, height = 6, seed = 123,
    color_by_annotation = COLOR_BY_ANNOTATION,
    annotations_file = file.path(results, "combined_node_scores.csv"),
    layout_by_group = LAYOUT_BY_GROUP, groups = NULL,
    group_layout = "connectivity") {
  if (!is.null(genes) && !is.null(gene_sets)) stop("Supply genes or gene_sets, not both")
  if (is.null(gene_sets)) gene_sets <- genes
  if (is.null(gene_sets)) stop("Supply genes or gene_sets")
  subset_network(results = results, gene_sets = gene_sets, out = out, title = title,
                 remove_isolates = remove_isolates, expand_one_hop = expand_one_hop,
                 color_node_groups = color_node_groups, node_groups = node_groups,
                 group_colors = group_colors, show_legend = show_legend,
                 font_family = font_family, width = width, height = height, seed = seed,
                 color_by_annotation = color_by_annotation, annotations_file = annotations_file,
                 layout_by_group = layout_by_group, groups = groups, group_layout = group_layout)
}
