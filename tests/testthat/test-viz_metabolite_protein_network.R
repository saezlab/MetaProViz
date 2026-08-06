library(testthat)
library(MetaProViz)

test_that("viz_metabolite_protein_network builds a network with metabolite names", {
    feature_metadata <- data.frame(
        Metabolite = c("Lactate", "Succinate"),
        hmdb = c("HMDB0000190", "1544"),
        stringsAsFactors = FALSE
    )

    metalinks_df <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0001544"),
        gene_symbol = c("HCAR1", "SLC13A3"),
        protein_type = c("receptor", "transporter"),
        protein_type_clean = c("Receptor", "Transporter"),
        transport_direction = c(NA, "in"),
        type = c("Ligand-Receptor", "Production-Degradation"),
        mode_of_regulation = c("Activating", NA),
        interaction_family = c("Receptor", "Transporter"),
        source = c("Example", "Example"),
        stringsAsFactors = FALSE
    )

    res <- viz_metabolite_protein_network(
        feature_metadata = feature_metadata,
        metalinks_df = metalinks_df,
        metabolite_col = "Metabolite",
        hmdb_col = "hmdb",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_true(is.list(res))
    expect_named(
        res,
        c("plot", "edges", "nodes", "matched_features", "unmatched_features", "saved_files")
    )
    expect_s3_class(res$plot, "ggplot")
    expect_gt(nrow(res$edges), 0)
    expect_gt(nrow(res$nodes), 0)
    expect_true(all(c("Lactate", "Succinate") %in% res$matched_features$metabolite))
    expect_length(res$saved_files, 0)
})

test_that("viz_metabolite_protein_network uses HMDB identifiers as fallback labels", {
    feature_metadata <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0000254; HMDB0000331"),
        stringsAsFactors = FALSE
    )

    metalinks_df <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0000254"),
        gene_symbol = c("HCAR1", "OXGR1"),
        type = c("Ligand-Receptor", "Ligand-Receptor"),
        mode_of_regulation = c("Activating", "Activating"),
        stringsAsFactors = FALSE
    )

    res <- viz_metabolite_protein_network(
        feature_metadata = feature_metadata,
        metalinks_df = metalinks_df,
        hmdb_col = "hmdb",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_true(all(c("HMDB0000190", "HMDB0000254") %in% res$matched_features$metabolite))
    expect_true(all(c("HMDB0000190", "HMDB0000254") %in% res$nodes$name))
    expect_true(any(res$edges$from == "HMDB0000190" | res$edges$to == "HMDB0000190"))
})

test_that("viz_metabolite_protein_network supports subsetted Metalinks tables", {
    feature_metadata <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0001544"),
        stringsAsFactors = FALSE
    )

    metalinks_df <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0001544"),
        gene_symbol = c("HCAR1", "SLC13A3"),
        type = c("Ligand-Receptor", "Production-Degradation"),
        stringsAsFactors = FALSE
    )

    res <- viz_metabolite_protein_network(
        feature_metadata = feature_metadata[1, , drop = FALSE],
        metalinks_df = metalinks_df[1, , drop = FALSE],
        hmdb_col = "hmdb",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_equal(nrow(res$matched_features), 1)
    expect_equal(sort(unique(c(res$edges$from, res$edges$to))), sort(c("HMDB0000190", "HCAR1")))
})

test_that("viz_metabolite_protein_network returns empty results when nothing matches", {
    feature_metadata <- data.frame(
        hmdb = c("HMDB0009999"),
        stringsAsFactors = FALSE
    )

    metalinks_df <- data.frame(
        hmdb = c("HMDB0000190"),
        gene_symbol = c("HCAR1"),
        type = c("Ligand-Receptor"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        res <- viz_metabolite_protein_network(
            feature_metadata = feature_metadata,
            metalinks_df = metalinks_df,
            hmdb_col = "hmdb",
            save_plot = NULL,
            print_plot = FALSE
        ),
        "No Metalinks interactions matched"
    )

    expect_null(res$plot)
    expect_equal(nrow(res$edges), 0)
    expect_equal(nrow(res$nodes), 0)
    expect_equal(nrow(res$matched_features), 1)
    expect_equal(nrow(res$unmatched_features), 1)
})

test_that("viz_metabolite_protein_network validates required columns", {
    feature_metadata <- data.frame(hmdb = "HMDB0000190", stringsAsFactors = FALSE)
    metalinks_df <- data.frame(hmdb = "HMDB0000190", stringsAsFactors = FALSE)

    expect_error(
        viz_metabolite_protein_network(
            feature_metadata = feature_metadata,
            metalinks_df = metalinks_df,
            hmdb_col = "hmdb",
            save_plot = NULL,
            print_plot = FALSE
        ),
        "`metalinks_df` must contain a `gene_symbol` column."
    )

    expect_error(
        viz_metabolite_protein_network(
            feature_metadata = feature_metadata,
            metalinks_df = data.frame(
                hmdb = "HMDB0000190",
                gene_symbol = "HCAR1",
                stringsAsFactors = FALSE
            ),
            metabolite_col = "Metabolite",
            hmdb_col = "hmdb",
            save_plot = NULL,
            print_plot = FALSE
        ),
        "`metabolite_col` not found"
    )
})
