library(testthat)
library(MetaProViz)

metalinks_info <- c(
    InputID = "hmdb",
    InputLabel = "Metabolite",
    PriorID = "hmdb",
    PriorTerm = "gene_symbol"
)

toy_metalinks <- function() {
    data.frame(
        hmdb = c("HMDB0000190", "HMDB0001544", "HMDB0000190"),
        gene_symbol = c("HCAR1", "SLC13A3", "SLC16A1"),
        protein_type_clean = c("gpcr", "transporter", "transporter"),
        interaction = c("Ligand-Receptor", "Transport_In", "Transport_Out"),
        direction = c("to_term", "to_metabolite", "to_term"),
        mode_of_regulation = c("Activating", "Binding", "Binding"),
        combined_score = c(900, 500, 700),
        stringsAsFactors = FALSE
    )
}


## viz_pk_network ------------------------------------------------------------

test_that("viz_pk_network builds a metabolite-term network", {
    feature_metadata <- data.frame(
        Metabolite = c("Lactate", "Succinate"),
        hmdb = c("HMDB0000190", "1544"),
        Log2FC = c(1.5, -2),
        stringsAsFactors = FALSE
    )

    res <- viz_pk_network(
        feature_metadata = feature_metadata,
        input_pk = toy_metalinks(),
        metadata_info = c(
            metalinks_info,
            MetaboliteColor = "Log2FC",
            TermColor = "protein_type_clean",
            EdgeColor = "interaction",
            EdgeLinetype = "mode_of_regulation",
            EdgeDirection = "direction"
        ),
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_named(res, c("DF", "Plot"))
    expect_named(res$DF, c("edges", "nodes", "matched_features", "unmatched_features"))
    expect_s3_class(res$Plot$pk_network, "ggplot")
    expect_equal(nrow(res$DF$edges), 3)
    expect_setequal(res$DF$edges$metabolite, c("Lactate", "Succinate"))
    expect_setequal(res$DF$edges$term, c("HCAR1", "SLC13A3", "SLC16A1"))
    expect_false(any(c(".directed", ".to_metabolite") %in% colnames(res$DF$edges)))

    nodes <- res$DF$nodes
    expect_equal(nodes$Log2FC[nodes$name == "Succinate"], -2)
    expect_true(is.na(nodes$Log2FC[nodes$name == "HCAR1"]))
    expect_equal(nodes$protein_type_clean[nodes$name == "HCAR1"], "gpcr")
    expect_equal(nodes$degree[nodes$name == "Lactate"], 2L)

    # The plot renders, including two fill scales via ggnewscale
    expect_silent(ggplot2::ggplot_build(res$Plot$pk_network))
})

test_that("viz_pk_network orients directed edges", {
    feature_metadata <- data.frame(
        Metabolite = c("Lactate", "Succinate"),
        hmdb = c("HMDB0000190", "HMDB0001544"),
        stringsAsFactors = FALSE
    )

    res <- viz_pk_network(
        feature_metadata = feature_metadata,
        input_pk = toy_metalinks(),
        metadata_info = c(metalinks_info, EdgeDirection = "direction"),
        save_plot = NULL,
        print_plot = FALSE
    )
    expect_s3_class(res$Plot$pk_network, "ggraph")

    edges <- MetaProViz:::.pk_network_edges(
        MetaProViz:::.pk_associations(
            feature_metadata,
            toy_metalinks(),
            c(metalinks_info, EdgeDirection = "direction"),
            id_type = "HMDB",
            id_sep = ";"
        )$associations,
        c(metalinks_info, EdgeDirection = "direction")
    )
    expect_true(all(edges$.directed))
    expect_equal(edges$.to_metabolite[edges$term == "SLC13A3"], TRUE)
    expect_equal(edges$.to_metabolite[edges$term == "HCAR1"], FALSE)
})

test_that("viz_pk_network works with exact IDs and a term_metadata table", {
    feature_metadata <- data.frame(
        Metabolite = c("Citrate", "Succinate", "Glucose"),
        KEGG = c("C00158", "C00042", "C00031"),
        stringsAsFactors = FALSE
    )
    kegg <- data.frame(
        term = c("TCA cycle", "TCA cycle", "Glycolysis", "Glycolysis"),
        MetaboliteID = c("C00158", "C00042", "C00031", "C00022"),
        stringsAsFactors = FALSE
    )
    ora <- data.frame(
        term = c("TCA cycle", "Glycolysis"),
        p.adjust = c(0.01, 0.2),
        Count = c(2, 1)
    )

    res <- viz_pk_network(
        feature_metadata = feature_metadata,
        input_pk = kegg,
        metadata_info = c(
            InputID = "KEGG",
            InputLabel = "Metabolite",
            PriorID = "MetaboliteID",
            PriorTerm = "term",
            TermColor = "p.adjust",
            TermSize = "Count"
        ),
        term_metadata = ora,
        id_type = "KEGG",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_equal(nrow(res$DF$edges), 3)
    nodes <- res$DF$nodes
    expect_equal(nodes$p.adjust[nodes$name == "TCA cycle"], 0.01)
    expect_equal(nodes$Count[nodes$name == "Glycolysis"], 1)
    expect_silent(ggplot2::ggplot_build(res$Plot$pk_network))
})

test_that("viz_pk_network splits multi-ID cells and labels by ID without labels", {
    feature_metadata <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0000254, HMDB0000331"),
        stringsAsFactors = FALSE
    )
    input_pk <- data.frame(
        hmdb = c("HMDB0000190", "HMDB0000254"),
        gene_symbol = c("HCAR1", "OXGR1"),
        stringsAsFactors = FALSE
    )

    res <- viz_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = c(InputID = "hmdb", PriorID = "hmdb", PriorTerm = "gene_symbol"),
        id_sep = ",",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_setequal(
        res$DF$nodes$name[res$DF$nodes$node_type == "Metabolite"],
        c("HMDB0000190", "HMDB0000254,HMDB0000331")
    )
    expect_equal(nrow(res$DF$matched_features), 2)
    expect_equal(res$DF$unmatched_features$id_normalized, "HMDB0000331")
})

test_that("metabolite and term nodes with the same name stay separate", {
    feature_metadata <- data.frame(Metabolite = "Choline", id = "1", stringsAsFactors = FALSE)
    input_pk <- data.frame(id = "1", term = "Choline", stringsAsFactors = FALSE)

    res <- viz_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = c(InputID = "id", InputLabel = "Metabolite", PriorID = "id", PriorTerm = "term"),
        id_type = "other",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_equal(nrow(res$DF$nodes), 2)
    expect_silent(ggplot2::ggplot_build(res$Plot$pk_network))
})

test_that("viz_pk_network drops mappings to columns with only missing values", {
    feature_metadata <- data.frame(Metabolite = "Lactate", hmdb = "HMDB0000190")
    input_pk <- toy_metalinks()
    input_pk$empty <- NA_character_

    expect_warning(
        res <- viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = input_pk,
            metadata_info = c(metalinks_info, EdgeColor = "empty"),
            save_plot = NULL,
            print_plot = FALSE
        ),
        "only missing values"
    )
    expect_s3_class(res$Plot$pk_network, "ggplot")
})

test_that("viz_pk_network returns empty results when nothing matches", {
    feature_metadata <- data.frame(hmdb = "HMDB0009999", stringsAsFactors = FALSE)

    expect_warning(
        res <- viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = toy_metalinks(),
            metadata_info = c(InputID = "hmdb", PriorID = "hmdb", PriorTerm = "gene_symbol"),
            save_plot = NULL,
            print_plot = FALSE
        ),
        "No feature"
    )

    expect_null(res$Plot$pk_network)
    expect_equal(nrow(res$DF$edges), 0)
    expect_equal(nrow(res$DF$matched_features), 0)
    expect_equal(nrow(res$DF$unmatched_features), 1)
})

test_that("viz_pk_network validates metadata_info", {
    feature_metadata <- data.frame(hmdb = "HMDB0000190", stringsAsFactors = FALSE)

    expect_error(
        viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = toy_metalinks(),
            metadata_info = c(InputID = "hmdb", PriorID = "hmdb"),
            save_plot = NULL,
            print_plot = FALSE
        ),
        "must contain: PriorTerm"
    )
    expect_error(
        viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = toy_metalinks(),
            metadata_info = metalinks_info,
            save_plot = NULL,
            print_plot = FALSE
        ),
        "`Metabolite`.*not found in `feature_metadata`"
    )
    expect_error(
        viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = toy_metalinks(),
            metadata_info = c(InputID = "hmdb", PriorID = "hmdb", PriorTerm = "gene_symbol", NodeColor = "x"),
            save_plot = NULL,
            print_plot = FALSE
        ),
        "Unknown `metadata_info` entries: NodeColor"
    )
    expect_error(
        viz_pk_network(
            feature_metadata = feature_metadata,
            input_pk = toy_metalinks(),
            metadata_info = c(InputID = "hmdb", PriorID = "hmdb", PriorTerm = "gene_symbol", EdgeWidth = "interaction"),
            save_plot = NULL,
            print_plot = FALSE
        ),
        "must name a numeric column"
    )
})


## viz_shared_pk_network -----------------------------------------------------

test_that("viz_shared_pk_network connects metabolites sharing terms", {
    feature_metadata <- data.frame(
        Metabolite = c("A", "B", "C", "D"),
        id = c("1", "2", "3", "4"),
        Log2FC = c(1, -1, 2, 0.5)
    )
    input_pk <- data.frame(
        id = c("1", "1", "1", "2", "2", "3", "4"),
        term = c("x", "y", "z", "x", "y", "z", "w")
    )
    info <- c(
        InputID = "id", InputLabel = "Metabolite", PriorID = "id", PriorTerm = "term",
        MetaboliteColor = "Log2FC", EdgeColor = "ignored_here"
    )

    res <- viz_shared_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = info,
        similarity = "jaccard",
        id_type = "other",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_named(res$DF, c("edges", "nodes", "associations", "matched_features", "unmatched_features"))
    edges <- res$DF$edges
    expect_equal(nrow(edges), 2)
    ab <- edges[edges$from == "A" & edges$to == "B", ]
    expect_equal(ab$shared, 2L)
    expect_equal(ab$jaccard, 2 / 3)
    expect_equal(ab$overlap_coefficient, 1)
    expect_equal(ab$weight, ab$jaccard)
    expect_equal(ab$shared_terms, "x; y")

    nodes <- res$DF$nodes
    expect_equal(nodes$n_terms[nodes$name == "A"], 3L)
    expect_equal(nodes$degree[nodes$name == "D"], 0L)
    expect_equal(nodes$Log2FC[nodes$name == "B"], -1)
    expect_silent(ggplot2::ggplot_build(res$Plot$shared_pk_network))

    res_threshold <- viz_shared_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = info[1:4],
        similarity = "jaccard",
        threshold = 0.5,
        id_type = "other",
        save_plot = NULL,
        print_plot = FALSE
    )
    expect_equal(nrow(res_threshold$DF$edges), 1)
})

test_that("viz_shared_pk_network plots metabolites without shared terms", {
    feature_metadata <- data.frame(Metabolite = c("A", "B"), id = c("1", "2"))
    input_pk <- data.frame(id = c("1", "2"), term = c("x", "y"))

    res <- viz_shared_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = c(InputID = "id", InputLabel = "Metabolite", PriorID = "id", PriorTerm = "term"),
        id_type = "other",
        save_plot = NULL,
        print_plot = FALSE
    )

    expect_equal(nrow(res$DF$edges), 0)
    expect_equal(nrow(res$DF$nodes), 2)
    expect_silent(ggplot2::ggplot_build(res$Plot$shared_pk_network))
})


## Helpers -------------------------------------------------------------------

test_that(".set_similarity matches the pairwise definitions", {
    sets <- list(a = c("1", "2", "3"), b = c("2", "3"), c = "4", d = character(0))

    shared <- MetaProViz:::.set_similarity(sets, "shared")
    expect_equal(shared["a", "b"], 2)
    expect_equal(shared["a", "c"], 0)

    jaccard <- MetaProViz:::.set_similarity(sets, "jaccard")
    expect_equal(jaccard["a", "b"], 2 / 3)
    expect_equal(jaccard["b", "a"], 2 / 3)
    expect_equal(jaccard["c", "d"], 0)
    expect_equal(unname(diag(jaccard)), rep(1, 4))

    overlap <- MetaProViz:::.set_similarity(sets, "overlap_coefficient")
    expect_equal(overlap["a", "b"], 1)
    expect_equal(overlap["a", "c"], 0)
})

test_that(".normalize_pk_ids normalises HMDB and trims other IDs", {
    expect_equal(
        MetaProViz:::.normalize_pk_ids(c("HMDB00123", "123", " hmdb0000123", "NA", "x"), "HMDB"),
        c("HMDB0000123", "HMDB0000123", "HMDB0000123", NA, NA)
    )
    expect_equal(
        MetaProViz:::.normalize_pk_ids(c(" C00031 ", "", "NA"), "KEGG"),
        c("C00031", NA, NA)
    )
})
