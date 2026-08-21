normalize_phosphosite_id <- function(x) {
  site_match <- stringr::str_match(
    stringr::str_trim(as.character(x)),
    "^(.+)_([STYsty])(\\d+)$"
  )

  dplyr::if_else(
    !is.na(site_match[, 1]),
    paste0(site_match[, 2], "_", stringr::str_to_upper(site_match[, 3]), site_match[, 4]),
    NA_character_
  )
}

phosphosite_parent <- function(x) {
  stringr::str_remove(as.character(x), "_[STY]\\d+$")
}

select_phosphosite_rows <- function(df, phosphosite_level = TRUE) {
  required_columns <- c("site_id", "logFC", "adj.P.Val")
  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0) {
    stop(
      "Phosphosite data are missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  selected <- df %>%
    dplyr::mutate(
      site_id = normalize_phosphosite_id(site_id),
      logFC = as.numeric(logFC),
      adj.P.Val = as.numeric(adj.P.Val),
      protein = phosphosite_parent(site_id),
      site_score = sign(logFC) * abs(logFC) * -log10(adj.P.Val + 1e-300)
    ) %>%
    dplyr::filter(
      !is.na(site_id), site_id != "",
      !is.na(logFC), !is.na(adj.P.Val)
    )

  if (isTRUE(phosphosite_level)) {
    selected <- selected %>%
      dplyr::group_by(site_id) %>%
      dplyr::slice_max(order_by = abs(site_score), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(node = site_id)
  } else {
    selected <- selected %>%
      dplyr::group_by(protein) %>%
      dplyr::slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(node = protein)
  }

  selected %>%
    dplyr::mutate(representative_site = site_id) %>%
    dplyr::relocate(node, representative_site, protein, site_id)
}

extract_signed_kinase_effects <- function(df) {
  kinase_column <- intersect(c("Kinase", "kinase", "hgnc_symbol"), names(df))
  effect_column <- intersect(
    c("effect_size", "kinase_effect", "signed_effect", "signed_activity", "activity", "z_score"),
    names(df)
  )

  if (length(kinase_column) == 0 || length(effect_column) == 0) {
    return(tibble::tibble(Symbol = character(), signed_effect = numeric()))
  }

  df %>%
    dplyr::transmute(
      Symbol = stringr::str_trim(as.character(.data[[kinase_column[[1]]]])),
      signed_effect = suppressWarnings(as.numeric(.data[[effect_column[[1]]]]))
    ) %>%
    dplyr::filter(!is.na(Symbol), Symbol != "") %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(
      signed_effect = if (all(is.na(signed_effect))) NA_real_ else mean(signed_effect, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_signed_kinase_effects <- function(kinase_data) {
  if (is.null(kinase_data) || nrow(kinase_data) == 0 ||
      !"signed_effect" %in% names(kinase_data)) {
    return(tibble::tibble(kinase = character(), kinase_effect = numeric()))
  }

  kinase_data %>%
    dplyr::transmute(
      kinase = as.character(name),
      signed_effect = as.numeric(signed_effect)
    ) %>%
    dplyr::group_by(kinase) %>%
    dplyr::summarise(
      kinase_effect = {
        known <- signed_effect[is.finite(signed_effect) & signed_effect != 0]
        if (length(known) == 0 || length(unique(sign(known))) != 1) NA_real_ else mean(known)
      },
      .groups = "drop"
    )
}

filter_kinase_phosphosite_sign <- function(
    edges,
    kinase_effects,
    phosphosite_effects,
    filter_sign = TRUE
) {
  if (!isTRUE(filter_sign) || is.null(edges) || nrow(edges) == 0) {
    return(edges)
  }

  required_edge_columns <- c("edge_type", "kinase", "phosphosite")
  missing_edge_columns <- setdiff(required_edge_columns, names(edges))
  if (length(missing_edge_columns) > 0) {
    stop(
      "Site-specific edges are missing required metadata column(s): ",
      paste(missing_edge_columns, collapse = ", ")
    )
  }

  kinase_effects <- kinase_effects %>%
    dplyr::select(kinase, kinase_effect) %>%
    dplyr::distinct(kinase, .keep_all = TRUE)

  phosphosite_effects <- phosphosite_effects %>%
    dplyr::transmute(
      phosphosite = as.character(site_id),
      phosphosite_effect = as.numeric(logFC)
    ) %>%
    dplyr::distinct(phosphosite, .keep_all = TRUE)

  edges %>%
    dplyr::left_join(kinase_effects, by = "kinase") %>%
    dplyr::left_join(phosphosite_effects, by = "phosphosite") %>%
    dplyr::mutate(
      .opposite_sign = edge_type == "kinase_to_phosphosite" &
        !is.na(kinase_effect) & is.finite(kinase_effect) & kinase_effect != 0 &
        !is.na(phosphosite_effect) & is.finite(phosphosite_effect) & phosphosite_effect != 0 &
        sign(kinase_effect) * sign(phosphosite_effect) < 0
    ) %>%
    dplyr::filter(!.opposite_sign) %>%
    dplyr::select(-kinase_effect, -phosphosite_effect, -.opposite_sign)
}

add_network_edge_metadata <- function(edges, source_name, edge_type = "protein_context") {
  edges %>%
    dplyr::mutate(
      source = .env$source_name,
      edge_type = .env$edge_type,
      kinase = NA_character_,
      phosphosite = NA_character_,
      parent_protein = NA_character_
    )
}

build_phosphoproteomic_network <- function(
    phosphosite_level,
    species,
    selected_phosphosite_rows,
    kinase_effects,
    include_string_context,
    filter_sign,
    generate_omnipath_site_edges_fn,
    generate_networkin_site_edges_fn,
    generate_phosphositeplus_site_edges_fn,
    generate_string_ppi_fn,
    generate_phuego_ppi_fn,
    generate_site_parent_edges_fn
) {
  if (isTRUE(phosphosite_level)) {
    observed_sites <- unique(selected_phosphosite_rows$site_id)

    phosphorylation_edges <- if (length(observed_sites) > 0) {
      dplyr::bind_rows(
        generate_omnipath_site_edges_fn(species, observed_sites),
        generate_networkin_site_edges_fn(species, observed_sites),
        generate_phosphositeplus_site_edges_fn(species, observed_sites)
      ) %>%
        filter_kinase_phosphosite_sign(
          kinase_effects = kinase_effects,
          phosphosite_effects = selected_phosphosite_rows,
          filter_sign = filter_sign
        )
    } else NULL

    string_context_edges <- if (isTRUE(include_string_context)) {
      generate_string_ppi_fn(species) %>%
        add_network_edge_metadata("STRING_context")
    } else NULL

    site_parent_edges <- if (isTRUE(include_string_context) && length(observed_sites) > 0) {
      generate_site_parent_edges_fn(observed_sites)
    } else NULL

    return(list(
      mode = "phosphosite_level",
      edges = dplyr::bind_rows(
        phosphorylation_edges,
        string_context_edges,
        site_parent_edges
      )
    ))
  }

  list(
    mode = "legacy_protein_level_STRING_phuEGO",
    edges = dplyr::bind_rows(
      generate_string_ppi_fn(species) %>%
        add_network_edge_metadata("STRING"),
      generate_phuego_ppi_fn(species) %>%
        add_network_edge_metadata("phuEGO")
    )
  )
}
