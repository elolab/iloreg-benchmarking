### Functions that are repeatedly called in the benchmarking ###

# ggplot2-based scatter plot visualization
PlotTwoDimensional <- function(two.dim.data = NULL,color.by=NULL,return.plot=FALSE,legend.title="group",xlab="dim1",ylab="dim2",title="",point.size=0.75)
{
    library(ggplot2)
    library(reshape2)
    library(scales)
    
    df <- as.data.frame(two.dim.data)
    if (!is.null(color.by)) {
        df$group <- color.by
        colnames(df) <- c("dim1","dim2","group")
        names(color.by) <- rownames(two.dim.data)
    } else {
        colnames(df) <- c("dim1","dim2")
        
    }
    if (class(color.by) == "numeric")
    {
        df$group[is.infinite(df$group)] <- NA
        df$type <- factor(is.na(df$group))
        p<-ggplot(df, aes(x=dim1, y=dim2)) +
            geom_point(size=point.size,aes(color=group)) +
            scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                   high = muted("blue"),name = legend.title,na.value = "gray20") +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
    } else if (is.null(color.by)) {
        p <- ggplot(df, aes(x=dim1, y=dim2)) + 
            geom_point(size=point.size) + 
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title)
    } else {
        
        two.dim.data_ <- two.dim.data
        rownames(two.dim.data_) <- names(color.by)
        cluster_centers <- lapply(levels(color.by),function(x) apply(two.dim.data_[names(color.by)[color.by==x],,drop=FALSE],2,median))
        cluster_centers <- do.call(rbind,cluster_centers)
        
        
        p<-ggplot(df, aes(x=dim1, y=dim2)) +
            geom_point(size=point.size,aes(color=group)) +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
            annotate("text", x = cluster_centers[,1], y = cluster_centers[,2], label = levels(color.by))
        
        
        
    }
    if (return.plot)
    {
        return(p)
    } else {
        print(p)
    }
}



# Wrapper for running ILoReg
RunILoReg <- function(data=NULL)
{
    # To install the correct version, run devtools::install_github("elolab/ILoReg",ref = "85196be6")
    library(ILoReg)
    
    iloreg_object <- CreateILoRegObject(normalized.data = data)
    set.seed(1)
    iloreg_object <- RunILoRegConsensus(iloreg_object,threads = 3,L = 200,k = 15,d = 0.3,r = 5,C = 0.3)
    set.seed(1) # To make the plots reproducable as well, set seed again because the first part in the manuscript's result was run in a cluster
    VisualizeConsensusInformation(iloreg_object,return.plot = F)
    iloreg_object <- ILoReg::RunPCA(iloreg_object,number.of.pcs = 50,scale = FALSE) # ILoReg:: to avoid overlap with Seurat
    PCAElbowPlot(iloreg_object,return.plot = F)
    iloreg_object <- ILoReg::RunTSNE(iloreg_object,perplexity=30) # ILoReg:: to avoid overlap with Seurat
    iloreg_object <- HierarchicalClustering(iloreg_object)
    iloreg_object <- CalculateSilhouetteInformation(iloreg_object,k.range = 2:50)
    SilhouetteCurve(iloreg_object)
    
    return(iloreg_object)
}


# Wrapper for running Seurat
RunSeurat <- function(raw.data=NULL, log.normalized.data=NULL, min.cells = 3, min.features = 200,nFeature_RNA_min = 200, nFeature_RNA_max = 2500, percent_mt = 5, nfeatures = 2000, dims = 1:15, resolution = 0.5,return.seurat.obj=FALSE)
{
    library(Seurat)
    
    if (class(raw.data) == "dgCMatrix")
    {
        raw.data <- as(raw.data, "dgTMatrix")
    }
    
    pbmc <- CreateSeuratObject(counts = raw.data, project = "pbmc3k", min.cells = min.cells, min.features = min.features)
    
    pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
    
    VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    CombinePlots(plots = list(plot1, plot2))
    
    print(nFeature_RNA_min)
    print(nFeature_RNA_max)
    print(percent_mt)
    
    
    
    pbmc <- subset(x = pbmc, subset = nFeature_RNA > get("nFeature_RNA_min"))
    pbmc <- subset(x = pbmc, subset = nFeature_RNA < get("nFeature_RNA_max"))
    pbmc <- subset(x = pbmc, subset = percent.mt < get("percent_mt"))
    
    pbmc <- NormalizeData(object = pbmc)
    
    if (!is.null(log.normalized.data))
    {
        rownames(log.normalized.data) <- gsub("_","-",rownames(log.normalized.data))
        cat("Replaced normalized data with log.normalized.data\n")
        pbmc[["RNA"]]@data <- log.normalized.data[rownames(pbmc[["RNA"]]@counts),colnames(pbmc[["RNA"]]@counts)]
    }
    
    pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = nfeatures)
    
    # Identify the 10 most highly variable genes
    top10 <- head(x = VariableFeatures(object = pbmc), 10)
    
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(object = pbmc)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
    
    
    all.genes <- rownames(x = pbmc)
    pbmc <- ScaleData(object = pbmc, features = all.genes)
    
    pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
    
    ElbowPlot(object = pbmc)
    
    pbmc <- FindNeighbors(object = pbmc, dims = dims)
    pbmc <- FindClusters(object = pbmc, resolution = resolution)
    
    pbmc <- Seurat::RunTSNE(object = pbmc, dims = dims)
    
    
    if (!return.seurat.obj)
    {
        return(list(clustering=Idents(object = pbmc),tsne=pbmc@reductions$tsne@cell.embeddings))
    } else {
        return(pbmc)
    }
    
}



# Wrapper for running SC3
RunSC3 <- function(raw.data=NULL,log.normalized.data=NULL,k=3,gene.filter=TRUE,n.cores=3)
{
    library(SingleCellExperiment)
    library(SC3)
    library(scater)
    
    ann <- as.data.frame(matrix(paste0("cell",1:ncol(raw.data)),ncol = 1,
                                dimnames = list(colnames(raw.data),"cell_type1")))
    
    sce <- SingleCellExperiment(
        assays = list(
            counts = as.matrix(raw.data),
            logcounts = as.matrix(log.normalized.data)
        ), 
        colData = ann
    )
    
    # define feature names in feature_symbol column
    rowData(sce)$feature_symbol <- rownames(sce)
    # remove features with duplicated names
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    
    
    sce <- sc3(sce, ks = k, biology = FALSE,gene_filter = gene.filter,n_cores=n.cores)
    
    return(sce)
    
}



# wrapper for running CIDR
RunCIDR <- function(raw.data=NULL,scale.factor=1e6,run.dimensionality.reduction=FALSE,return.cidr.object=TRUE)
{
    library(cidr)
    library(Matrix)
    
    
    raw.data <- as.matrix(raw.data)
    
    # This runs a regular prcomp analysis, as in the tutorial
    if (run.dimensionality.reduction)
    {
        
        priorTPM <- 1
        pan10 <- raw.data[rowSums(raw.data)>10,]
        pan10_lcpm <- log2(t(t(pan10)/colSums(pan10))*scale.factor+priorTPM)
        PC_lcpm <- prcomp(t(pan10_lcpm))
        
    }
    
    sData <- scDataConstructor(raw.data)
    
    sData <- determineDropoutCandidates(sData)
    
    sData <- wThreshold(sData)
    
    sData <- scDissim(sData)
    
    sData <- scPCA(sData)
    
    sData <- nPC(sData)
    
    sData <- scCluster(sData)
    

    if (!return.cidr.object)
    {
        
        if (run.dimensionality.reduction)
        {
            return(list(pca=sData@PC[,c(1,2)],
                        clustering=factor(sData@clusters)))
        } else {
            return(factor(sData@clusters))
        }
        
    } else {
        if (run.dimensionality.reduction)
        {
            return(list(pca=sData@PC[,c(1,2)],
                        cidr_out=sData))
        } else {
            return(sData)
        }
        
    }
}



# Wrapper for running RaceID3
RunRaceID <- function(raw.data=NULL,normalized.data=NULL,max.clusters=30)
{
    library(RaceID)
    
    sc <- SCseq(raw.data)
    
    sc <- filterdata(sc,mintotal=1)
    
    # if we want to force the tool to use specific normalized data
    # not used in this benchmarking, though
    if(!is.null(normalized.data))
    {
        sc@ndata <- normalized.data
        sc@genes <- rownames(normalized.data)
    }
    
    sc <- compdist(sc,metric="pearson")
    
    sc <- clustexp(sc,clustnr = max.clusters)
    
    plotsaturation(sc,disp=FALSE)
    
    plotsaturation(sc,disp=TRUE)
    
    plotjaccard(sc)
    
    sc <- findoutliers(sc)
    
    plotbackground(sc)
    
    plotsensitivity(sc)
    
    plotoutlierprobs(sc)
    
    clustheatmap(sc)
    
    sc <- comptsne(sc)
    
    
    return(sc)
    
}
