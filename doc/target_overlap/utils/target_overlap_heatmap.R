# heatmap plotting function for enrichments and depletions of target overlaps

require(pheatmap)
require(RColorBrewer)
require(dplyr)

plot.enrichments <- function(df,
                             categories,
                             exprThresh = 10,
                             cluster.rows = FALSE,
                             cluster.cols = FALSE,
                             row.order = NULL,
                             column.order = NULL,
                             title = "Target Enrichments") {
    # a "well-formatted" data frame, specific to the target_overlap enrichments created in
    # target_overlap_cell_type_specific_fisher.R

    # filter for only "category" mirnas
    df <- dplyr::filter(df, category %in% categories)

    # filter for high expressers: expression in gzcp takes precidence for expression thresholds
    if ("NeuroGZ" %in% categories | "NeuroCP" %in% categories) {
        df <- dplyr::filter(df, baseMean.gzcp > exprThresh)
    } else {
        df <- dplyr::filter(df, baseMean.gw > exprThresh)
    }

    # remove mirnas with no enrichments or depletions
    df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0 |
                            rowSums(dplyr::select(df, starts_with("FTdepleted"))) > 0)

    # set l2fold change to 0 if not significant
    df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)
    df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

    # convert to non-tibble data frame
    df <- as.data.frame(df)
    rownames(df) <- df$Name

    # matrix of odds ratios for plotting
    mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
    # significant enrichments only
    mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))
    mat.depleted <- as.matrix(dplyr::select(df, starts_with("FTdepleted")))

    # relabel column names
    colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
    colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
    colnames(mat.depleted) <- sapply(strsplit(colnames(mat.depleted), "_"), `[`, 2)

    # check column names are the same across the three matrices
    if ( all(colnames(mat.odds) != colnames(mat.enriched)) |
         all(rownames(mat.odds) != rownames(mat.enriched)) |
         all(colnames(mat.odds) != colnames(mat.depleted)) |
         all(rownames(mat.odds) != rownames(mat.depleted)) ) {
        stop("Matrix column names not the same.")
    }

    # df for cluster labels
    CellType <- factor(c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron",
                         "InNeuron", "Progenitor", "BulkDE", "BulkDE", "Other", "BulkDE", "BulkDE", "Other",
                         "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor"),
                       levels = c("BulkDE", "Progenitor", "ExNeuron", "InNeuron", "Other"), ordered = TRUE)
    df.cluster <- data.frame(CellType = CellType)
    #df.cluster$CellType <- as.factor(df.cluster$CellType)
    rownames(df.cluster) <- colnames(mat.odds)

    # plotting matrix
    mat <- mat.odds
    # if not enriched, set odds to 1
    mat[!mat.enriched & !mat.depleted] <- 1
    mat[is.na(mat)] <- 1
    mat[is.infinite(mat)] <- 1
    mat[mat == 0] <- 1

    # transform odds ratios
    mat <- log10(mat)

    # order matrix
    if (is.null(row.order) &
        ("NeuroGZ" %in% categories | "NeuroCP" %in% categories)) {
        row.order <- order(df$log2FoldChange.gzcp)
    } else {
        row.order <- order(df$log2FoldChange.gw)
    }
    if (is.null(column.order)) {
        column.order <- rownames(df.cluster)[order(df.cluster$CellType)]
    }
    mat <- mat[row.order, column.order]

    # df for row labels
    df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                          L2FC_GW = df$log2FoldChange.gw,
                          L2Expression = log2(df$baseMean.gzcp),
                          Source = df$source)
    rownames(df.mirs) <- df$Name

    #cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
    #cols <- colorRampPalette(c("white", "red3"))(10)

    reds <- colorRampPalette(c("white", "red3"))(15)
    blues <- colorRampPalette(c("white", "blue3"))(15)

    cols <- c(rev(blues), reds)

    ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                       L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                       L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                       CellType = c(Other = "#1B9E77",
                                    ExNeuron = "#D95F02",
                                    InNeuron = "#7570B3",
                                    Progenitor = "#E7298A",
                                    BulkDE = "grey50"),
                       Source = c(miRBase_v22 = "#66A61E",
                                  Friedlander2014 = "#E6AB02",
                                  miRDeep2 = "#BF5B17",
                                  miRge = "#386CB0",
                                  Nowakowski2018 = "grey70"))

    #rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

    # create heatmap
    # row cluster order: bulkDE, progenitors, exNeuron, inneuron, other
    p <- pheatmap(mat = mat,
                  color = cols,
                  cluster_rows = cluster.rows,
                  cluster_cols = cluster.cols,
                  show_colnames = TRUE,
                  show_rownames = TRUE,
                  annotation_col = df.cluster,
                  annotation_row = df.mirs,
                  na_col = "grey50",
                  breaks = seq(from=-1.5, to=1.5, by=0.1),
                  annotation_colors = ann_colors,
                  main = paste(title, "(>", exprThresh, "mean counts)"),
                  fontsize_row = 4,
                  gaps_col = c(4,9,14,16))

    print(p)
}

plot.lg.enrichments <- function(df,
                                categories,
                                exprThresh = 10,
                                cluster.rows = FALSE,
                                cluster.cols = FALSE,
                                row.order = NULL,
                                column.order = NULL,
                                title = "Target Enrichments") {

    # filter for only "category" mirnas
    df <- dplyr::filter(df, category %in% categories)

    # filter for high expressers: expression in gzcp takes precidence for expression thresholds
    if ("NeuroGZ" %in% categories | "NeuroCP" %in% categories) {
        df <- dplyr::filter(df, baseMean.gzcp > exprThresh)
    } else {
        df <- dplyr::filter(df, baseMean.gw > exprThresh)
    }

    # remove mirnas with no enrichments or depletions
    df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("LGenriched"))) > 0 |
                            rowSums(dplyr::select(df, starts_with("LGdepleted"))) > 0)

    # set l2fold change to 0 if not significant
    df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)
    df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

    # convert to non-tibble data frame
    df <- as.data.frame(df)
    rownames(df) <- df$Name

    # matrix of odds ratios for plotting
    mat.odds <- as.matrix(dplyr::select(df, starts_with("LGoddsRatio")))
    # significant enrichments only
    mat.enriched <- as.matrix(dplyr::select(df, starts_with("LGenriched")))
    mat.depleted <- as.matrix(dplyr::select(df, starts_with("LGdepleted")))

    # relabel column names
    colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
    colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
    colnames(mat.depleted) <- sapply(strsplit(colnames(mat.depleted), "_"), `[`, 2)

    # check column names are the same across the three matrices
    if ( all(colnames(mat.odds) != colnames(mat.enriched)) |
         all(rownames(mat.odds) != rownames(mat.enriched)) |
         all(colnames(mat.odds) != colnames(mat.depleted)) |
         all(rownames(mat.odds) != rownames(mat.depleted)) ) {
        stop("Matrix column names not the same.")
    }

    # df for cluster labels
    CellType <- factor(c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron",
                         "InNeuron", "Progenitor", "BulkDE", "BulkDE", "Other", "BulkDE", "BulkDE", "Other",
                         "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor"),
                       levels = c("BulkDE", "Progenitor", "ExNeuron", "InNeuron", "Other"), ordered = TRUE)
    df.cluster <- data.frame(CellType = CellType)
    #df.cluster$CellType <- as.factor(df.cluster$CellType)
    rownames(df.cluster) <- colnames(mat.odds)

    # plotting matrix
    mat <- mat.odds
    # if not enriched, set odds to 1
    mat[!mat.enriched & !mat.depleted] <- 1
    mat[is.na(mat)] <- 1
    mat[is.infinite(mat)] <- 1
    mat[mat == 0] <- 1

    # transform odds ratios
    mat <- log10(mat)

    # order matrix
    if (is.null(row.order) &
        ("NeuroGZ" %in% categories | "NeuroCP" %in% categories)) {
        row.order <- order(df$log2FoldChange.gzcp)
    } else {
        row.order <- order(df$log2FoldChange.gw)
    }
    if (is.null(column.order)) {
        column.order <- rownames(df.cluster)[order(df.cluster$CellType)]
    }
    mat <- mat[row.order, column.order]

    # df for row labels
    df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                          L2FC_GW = df$log2FoldChange.gw,
                          L2Expression = log2(df$baseMean.gzcp),
                          Source = df$source)
    rownames(df.mirs) <- df$Name

    #cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
    #cols <- colorRampPalette(c("white", "red3"))(10)

    reds <- colorRampPalette(c("white", "red3"))(15)
    blues <- colorRampPalette(c("white", "blue3"))(15)

    cols <- c(rev(blues), reds)

    ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                       L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                       L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                       CellType = c(Other = "#1B9E77",
                                    ExNeuron = "#D95F02",
                                    InNeuron = "#7570B3",
                                    Progenitor = "#E7298A",
                                    BulkDE = "grey50"),
                       Source = c(miRBase_v22 = "#66A61E",
                                  Friedlander2014 = "#E6AB02",
                                  miRDeep2 = "#BF5B17",
                                  miRge = "#386CB0",
                                  Nowakowski2018 = "grey70"))

    #rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

    # create heatmap
    # row cluster order: bulkDE, progenitors, exNeuron, inneuron, other
    p <- pheatmap(mat = mat,
                  color = cols,
                  cluster_rows = cluster.rows,
                  cluster_cols = cluster.cols,
                  show_colnames = TRUE,
                  show_rownames = TRUE,
                  annotation_col = df.cluster,
                  annotation_row = df.mirs,
                  na_col = "grey50",
                  breaks = seq(from=-1.5, to=1.5, by=0.1),
                  annotation_colors = ann_colors,
                  main = paste(title, "(>", exprThresh, "mean counts)"),
                  fontsize_row = 4,
                  gaps_col = c(4,9,14,16))

    print(p)
}
