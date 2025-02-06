read10x_mtx <- function(run, suffix, min_counts=1) {
    #min_counts initially 100
    barcode.loc <- list.files(run, pattern = 'barcodes.tsv(.gz)?', full.names = TRUE)
    gene.loc <- list.files(run, pattern = 'features.tsv(.gz)?', full.names = TRUE)
    matrix.loc <- list.files(run, pattern = 'matrix.mtx(.gz)?', full.names = TRUE)
    

    data <- readMM(file = matrix.loc) %>% as("dgCMatrix")# %>% Matrix::t()
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)    
    
    if (!missing(suffix)) {
        cell.names <- paste(cell.names, suffix, sep = "_")
    }
    gene.names <- fread(gene.loc, header = FALSE)$V2
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    data <- as(data, "dgCMatrix")
    data <- data[, Matrix::colSums(data) >= min_counts]
    data <- data[which(!is.na(row.names(data))), ]
    data <- as(sumOverRowNames(data), "dgCMatrix")
    return(data)
}

sumOverRowNames <- function(X) {
    name_factors <- factor(row.names(X))
    res <- presto::sumGroups(X, name_factors)
    row.names(res) <- levels(name_factors)#[1:nrow(res)]
    colnames(res) <- colnames(X)
    return(res)
}

bin_95 <- function (data_df, xvar, yvars, num.bin, .include.lowest = TRUE) 
{
    data_df <- data.frame(data_df)
    x <- data_df[[xvar]]
    .breaks <- unique(quantile(x, probs = seq(0, 1, length.out = num.bin)))
    bin_means <- 0.5 * (head(.breaks, -1) + tail(.breaks, -1))
    bins_freq <- cut(x, .breaks, include.lowest = .include.lowest)
    levels(bins_freq) <- bin_means
    xmeans <- lapply(split(x, as.integer(bins_freq)), mean) %>% 
        as.numeric
    res_df <- Reduce(rbind, lapply(yvars, function(yvar) {
        y <- data_df[[yvar]]
        ymeans <- lapply(split(y, as.integer(bins_freq)), mean) %>% 
            as.numeric
        ysd <- lapply(split(y, as.integer(bins_freq)), sd) %>% 
            as.numeric
        data.frame(xval = xmeans) %>% cbind(Reduce(rbind, lapply(split(y, 
            as.integer(bins_freq)), function(.x) quantile(.x, 
            c(0.05, 0.95)))) %>% data.frame() %>% dplyr::mutate(symbol = yvar)) %>% 
            cbind(yval = ymeans, ysd = ysd)
    })) %>% data.frame()
    return(res_df)
}                                         


scDblFinder_par <- function(counts, library_ids, logcounts=NULL, ncore=1) {
    if (is.null(logcounts)) {
        logcounts <- singlecellmethods::normalizeData(counts, 1e4, 'log')
    }
    logcounts_list <- split(seq_len(ncol(counts)), library_ids) %>% 
        map(function(idx) {
            logcounts[, idx]
        })

    ## TODO: pick variable genes 
    
    counts_list <- split(seq_len(ncol(counts)), library_ids) %>% 
        map(function(idx) {
            counts[, idx]
        })

    if (ncore == 1) {
        future::plan(sequential)
    } else if (ncore %in% c(0, Inf)) {
        ncore <- availableCores()
        future::plan(multiprocess)
    } else {
        ## TODO: remove this assignment that pollutes the global environment 
        .ncore <<- ncore
        future::plan(future::multiprocess(workers = .ncore))
    }
    
    future_map2(logcounts_list, counts_list, function(.logcounts, .counts) {
            sce <- SingleCellExperiment(list(counts = .counts,logcounts = .logcounts))
            sce <- scDblFinder(sce, verbose=FALSE)
            as.data.frame(sce@colData) %>% 
                tibble::rownames_to_column('CellID')
        }) %>% 
        bind_rows(.id = 'LibraryID')
}
