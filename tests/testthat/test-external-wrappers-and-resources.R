test_that("config and logging wrappers delegate the MetaProViz package name", {
  calls <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    save_config = function(path, title, local, pkg) {
      calls$save <- list(path = path, title = title, local = local, pkg = pkg)
      NULL
    },
    load_config = function(path, title, user, pkg, ...) {
      calls$load <- list(path = path, title = title, user = user, pkg = pkg)
      list(ok = TRUE)
    },
    reset_config = function(save, reset_all, pkg) {
      calls$reset <- list(save = save, reset_all = reset_all, pkg = pkg)
      list(ok = TRUE)
    },
    config_path = function(user, pkg) {
      calls$config_path <- list(user = user, pkg = pkg)
      "mock-config.yml"
    },
    logfile = function(pkg) {
      calls$logfile <- pkg
      "mock.log"
    },
    read_log = function(pkg) {
      calls$read_log <- pkg
      "mock log"
    },
    set_loglevel = function(level, target, pkg) {
      calls$loglevel <- list(level = level, target = target, pkg = pkg)
      NULL
    },
    .package = "MetaProViz"
  )

  expect_null(metaproviz_save_config(path = "x.yml", title = "unit", local = TRUE))
  expect_identical(calls$save$pkg, "MetaProViz")

  expect_equal(metaproviz_load_config(path = "x.yml"), list(ok = TRUE))
  expect_identical(calls$load$pkg, "MetaProViz")

  expect_equal(metaproviz_reset_config(save = TRUE), list(ok = TRUE))
  expect_identical(calls$reset$pkg, "MetaProViz")

  expect_identical(metaproviz_config_path(user = TRUE), "mock-config.yml")
  expect_identical(calls$config_path$pkg, "MetaProViz")

  expect_identical(metaproviz_logfile(), "mock.log")
  expect_identical(calls$logfile, "MetaProViz")

  expect_identical(metaproviz_log(), "mock log")
  expect_identical(calls$read_log, "MetaProViz")

  expect_null(metaproviz_set_loglevel("trace", target = "console"))
  expect_identical(calls$loglevel$pkg, "MetaProViz")
})

test_that("metsigdb_kegg, reactome, wikipathways and macdb reshape mocked resource data", {
  testthat::local_mocked_bindings(
    kegg_link = function(...) data.frame(id_a = "path:hsa00010", id_b = "cpd:C00031", stringsAsFactors = FALSE),
    kegg_list = function(type) {
      if (identical(type, "pathway")) {
        data.frame(id = "hsa00010", name = "Glycolysis", stringsAsFactors = FALSE)
      } else {
        data.frame(id = "C00031", name = "Glucose; D-Glucose", stringsAsFactors = FALSE)
      }
    },
    kegg_conv = function(...) data.frame(id_a = "123", id_b = "C00031", stringsAsFactors = FALSE),
    reactome_chebi = function(...) data.frame(pathway_id = "R-HSA-1", pathway_name = "Reactome 1", chebi_ids = "CHEBI:1, CHEBI:2", stringsAsFactors = FALSE),
    wikipathways_metabolites_sparql = function(...) {
      data.frame(
        pathway_id = "WP1",
        pathway_name = "WP name",
        pathway_url = "http://example.com",
        n_metabolites_in_pathway = 2,
        metabolites = "CHEBI:1; CHEBI:2",
        stringsAsFactors = FALSE
      )
    },
    macdb_metabolite_cancer_associations = function(...) {
      data.frame(
        metabolite_pubchem_cid = c("123", "123"),
        metabolite_name = c("Lactate", "Lactate"),
        study_cancer_type = c("Cancer", "Cancer"),
        study_cancer_subtype = c("Subtype", "Subtype"),
        study_pubmed_id = c("1", "2"),
        association_pvalue = c(0.01, 0.2),
        association_log2fc = c(1, 2),
        Cohort_id = c("C1", "C2"),
        study_trait_onto_id = c("T1", "T1"),
        study_cancer_doid = c("D1", "D1"),
        study_tissue = c("Liver", "Liver"),
        association_case_concentration = c(1, 2),
        association_control_concentration = c(0.5, 0.4),
        has_pValue = c(TRUE, TRUE),
        has_log2FC = c(TRUE, TRUE),
        has_caseConcentration = c(TRUE, TRUE),
        has_controlConcentration = c(TRUE, TRUE),
        has_bothConcentrations = c(TRUE, TRUE),
        stringsAsFactors = FALSE
      )
    },
    .package = "MetaProViz"
  )

  kegg <- metsigdb_kegg(exclude_metabolites = NULL)
  expect_true(all(c("term", "Metabolite", "MetaboliteID", "Description") %in% colnames(kegg)))

  reactome <- metsigdb_reactome(exclude_metabolites = NULL)
  expect_true("chebi_id" %in% colnames(reactome))

  wp <- metsigdb_wikipathways(exclude_metabolites = NULL)
  expect_true("metabolite_id" %in% colnames(wp))

  macdb <- metsigdb_macdb(exclude_metabolites = NULL)
  expect_true(all(c("term", "Metabolite_PubchemID", "association_score") %in% colnames(macdb)))
})

test_that("metsigdb_chemicalclass and metsigdb_metalinks work with mocked backends", {
  testthat::local_mocked_bindings(
    user_cache_dir = function(...) tempfile("mpv-cache"),
    ramp_table = function(table, version = NULL) {
      if (identical(table, "metabolite_class")) {
        data.frame(
          ramp_id = c("R1", "R1"),
          class_source_id = c("hmdb:HMDB0000001", "hmdb:HMDB0000001"),
          class_level_name = c("ClassyFire_class", "ClassyFire_super_class"),
          class_name = c("Organic acids", "Lipids"),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          ramp_id = "R1",
          a = "x",
          b = "y",
          c = "z",
          d = 1,
          e = 2,
          f = 3,
          g = 4,
          h = 5,
          common_name = "Lactate",
          stringsAsFactors = FALSE
        )
      }
    },
    metalinksdb_sqlite = function(...) structure(list(), class = "mock_con"),
    dbListTables = function(con) names(mpv_minimal_metalinks_tables()),
    dbDisconnect = function(con) TRUE,
    dbGetQuery = function(con, query) {
      tables <- mpv_minimal_metalinks_tables()
      for (nm in names(tables)) {
        if (grepl(nm, query, fixed = TRUE)) {
          return(tables[[nm]])
        }
      }
      stop("Unexpected query: ", query)
    },
    .package = "MetaProViz"
  )

  chemical_class <- metsigdb_chemicalclass(save_table = NULL, exclude_metabolites = NULL)
  expect_true(all(c("class_source_id", "common_name") %in% colnames(chemical_class)))

  metalinks <- metsigdb_metalinks(save_table = NULL, exclude_metabolites = NULL)
  expect_true(all(c("hmdb", "gene_symbol", "type", "interaction_family") %in% colnames(metalinks)))
})
