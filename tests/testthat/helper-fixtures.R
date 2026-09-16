mpv_load_dataset <- function(name) {
  data(list = name, package = "MetaProViz", envir = environment())
  get(name, envir = environment())
}

mpv_intracell_fixture <- function(n = 12) {
  raw <- mpv_load_dataset("intracell_raw")
  raw <- raw[seq_len(min(n, nrow(raw))), , drop = FALSE]
  df <- tibble::column_to_rownames(raw, "Code")

  list(
    data = df[, -c(1:3), drop = FALSE],
    metadata_sample = df[, c(1:3), drop = FALSE],
    metadata_info = c(
      Conditions = "Conditions",
      Biological_Replicates = "Biological_Replicates",
      Analytical_Replicates = "Analytical_Replicates"
    )
  )
}

mpv_intracell_se_fixture <- function() {
  mpv_load_dataset("intracell_raw_se")
}

mpv_intracell_dma_fixture <- function(n = 20) {
  raw <- mpv_load_dataset("intracell_dma")
  raw[seq_len(min(n, nrow(raw))), , drop = FALSE]
}

mpv_core_dma_fixture <- function(n = 20) {
  intra <- mpv_intracell_dma_fixture(n)

  data.frame(
    Metabolite = intra$Metabolite,
    `Log2(Distance)` = seq(-1.5, 1.5, length.out = nrow(intra)),
    p.adj = seq(0.001, 0.05, length.out = nrow(intra)),
    core = rep(c("Consumption", "Released"), length.out = nrow(intra)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


mpv_medium_fixture <- function(n = 20) {
  raw <- mpv_load_dataset("medium_raw")
  raw <- raw[seq_len(min(n, nrow(raw))), , drop = FALSE]
  df <- tibble::column_to_rownames(raw, "Code")

  list(
    data = df[, -c(1:3), drop = FALSE],
    metadata_sample = df[, c(1:3), drop = FALSE],
    metadata_info = c(
      Conditions = "Conditions",
      Biological_Replicates = "Biological_Replicates",
      core_norm_factor = "GrowthFactor",
      core_media = "blank"
    )
  )
}

mpv_simple_pathway_fixture <- function() {
  data.frame(
    term = c("Pathway A", "Pathway A", "Pathway B", "Pathway B"),
    Metabolite = c("M1", "M2", "M2", "M3"),
    gene = c("M1", "M2", "M2", "M3"),
    Description = c("Desc A", "Desc A", "Desc B", "Desc B"),
    stringsAsFactors = FALSE
  )
}

mpv_cluster_input_fixture <- function() {
  data.frame(
    cluster = c("cluster1", "cluster1", "cluster2"),
    BG_method = c("TRUE", "TRUE", "TRUE"),
    stringsAsFactors = FALSE,
    row.names = c("M1", "M2", "M3")
  )
}

mpv_standard_ora_input_fixture <- function() {
  data.frame(
    p.adj = c(0.001, 0.002, 0.2),
    t.val = c(12, -11, 0.5),
    stringsAsFactors = FALSE,
    row.names = c("M1", "M2", "M3")
  )
}

mpv_seed_edge_table <- function() {
  tibble::tibble(
    id1 = c("HMDB0000001", "C00001", "C00001", "1"),
    type1 = c("HMDB", "KEGG", "KEGG", "CHEBI"),
    id2 = c("C00001", "HMDB0000001", "1", "C00001"),
    type2 = c("KEGG", "HMDB", "CHEBI", "KEGG")
  )
}

mpv_seed_input_data <- function() {
  data.frame(
    metabolite = c("full", "partial", "complete_hmdb", "complete_chebi"),
    HMDB = c("HMDB1", "HMDB2; HMDB3", "HMDB4; HMDB5", NA),
    KEGG = c("C00001", "C00002", NA, "C00003"),
    CHEBI = c("CHEBI:1", NA, NA, "CHEBI:3; CHEBI:4"),
    PUBCHEM = NA_character_,
    stringsAsFactors = FALSE
  )
}

mpv_compare_pk_fixture <- function() {
  list(
    KEGG = data.frame(term = c("T1", "T1", "T2"), MetaboliteID = c("HMDB0000001", "HMDB0000002", "HMDB0000003")),
    Reactome = data.frame(term = c("R1", "R2"), MetaboliteID = c("HMDB0000002", "HMDB0000004"))
  )
}

mpv_translation_fixture <- function() {
  data.frame(
    MetaboliteID = c("A", "B", "C"),
    term = c("P1", "P1", "P2"),
    stringsAsFactors = FALSE
  )
}

mpv_mock_translate_ids <- function(data, ...) {
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  n <- nrow(out)
  if (!"hmdb" %in% colnames(out)) out$hmdb <- paste0("HMDB", seq_len(n))
  if (!"kegg" %in% colnames(out)) out$kegg <- paste0("C", sprintf("%05d", seq_len(n)))
  if (!"chebi" %in% colnames(out)) out$chebi <- paste0("CHEBI:", seq_len(n))
  if (!"pubchem" %in% colnames(out)) out$pubchem <- as.character(seq_len(n))
  out
}

mpv_minimal_metalinks_tables <- function() {
  list(
    interactions = data.frame(
      hmdb = c("HMDB0000190", "HMDB0001544", "HMDB0000002"),
      uniprot = c("P12345", "P23456", "P34567"),
      type = c("lr", "lr", "lr"),
      mor = c(1, 1, 1),
      transport_direction = c(NA_character_, NA_character_, NA_character_),
      protein_type = c("gpcr", "transporter", "enzyme"),
      source = c("CellPhoneDB", "CellPhoneDB", "CellPhoneDB"),
      experiment_score = c(1, 1, 1),
      combined_score = c(900, 900, 900),
      stringsAsFactors = FALSE
    ),
    cell_location = data.frame(hmdb = "HMDB0000190", cell_location = "Extracellular", stringsAsFactors = FALSE),
    tissue_location = data.frame(hmdb = "HMDB0000190", tissue_location = "Liver", stringsAsFactors = FALSE),
    biospecimen_location = data.frame(hmdb = "HMDB0000190", biospecimen_location = "Blood", stringsAsFactors = FALSE),
    disease = data.frame(hmdb = "HMDB0000190", disease = "Cancer", stringsAsFactors = FALSE),
    pathway = data.frame(hmdb = "HMDB0000190", pathway = "Glycolysis", stringsAsFactors = FALSE),
    metabolites = data.frame(
      hmdb = c("HMDB0000190", "HMDB0001544", "HMDB0000002"),
      metabolite = c("Lactate", "Succinate", "Acetate"),
      stringsAsFactors = FALSE
    ),
    proteins = data.frame(
      uniprot = c("P12345", "P23456", "P34567"),
      gene_symbol = c("HCAR1", "SLC13A3", "ACSS2"),
      protein_name = c("Hydroxycarboxylic acid receptor 1", "Solute carrier family 13 member 3", "Acyl-CoA synthetase short-chain family member 2"),
      stringsAsFactors = FALSE
    )
  )
}
