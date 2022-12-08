# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(biomaRt)
library(hgu133plus2.db)
library(lumiHumanAll.db)

# load series and platform data from GEO

gset <- getGEO("GSE45827", GSEMatrix = TRUE, AnnotGPL=TRUE)
gene_ids <- gset$GSE45827_series_matrix.txt.gz@featureData$ID

if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11212211111222212112221211111111221112221122121122",
               "211211221111111112222XXXXXXXXXXXXXX000000000003344",
               "34434433443344433443334433444344334433333344344334",
               "33444")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# filter genes only with genes of interest
mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "ensembl_gene_id",
    "affy_hg_u133_plus_2",
    "hgnc_symbol", "description","refseq_mrna", "refseq_ncrna"),
  filter = "affy_hg_u133_plus_2",
  values = gene_ids,
  uniqueRows=TRUE)
gene.set.of.interest <- c("SNORA1","SNORA12","SNORA14B","SNORA16A","SNORA21","SNORA23","SNORA24","SNORA32","SNORA38","SNORA44","SNORA48","SNORA49","SNORA52","SNORA53","SNORA57","SNORA61","SNORA63","SNORA64","SNORA65","SNORA70","SNORA71C","SNORA73A","SNORA73B","SNORA75","SNORA78","SNORA8","SNORD10","SNORD102","SNORD104","SNORD105","SNORD105B","SNORD108","SNORD110","SNORD111B","SNORD113-3","SNORD113-4","SNORD114-1","SNORD114-13","SNORD114-14","SNORD114-19","SNORD114-20","SNORD114-21","SNORD115-23","SNORD115-32","SNORD116-13","SNORD119","SNORD12","SNORD12B","SNORD13","SNORD14A","SNORD14C","SNORD14D","SNORD15B","SNORD16","SNORD17","SNORD18A","SNORD1B","SNORD20","SNORD22","SNORD26","SNORD28","SNORD29","SNORD32A","SNORD33","SNORD34","SNORD35A","SNORD36A","SNORD36B","SNORD38A","SNORD3A","SNORD3D","SNORD41","SNORD42A","SNORD42B","SNORD44","SNORD45A","SNORD47","SNORD49B","SNORD4B","SNORD5","SNORD50A","SNORD52","SNORD53","SNORD54","SNORD55","SNORD56","SNORD57","SNORD58A","SNORD58C","SNORD59A","SNORD60","SNORD63","SNORD64","SNORD65","SNORD68","SNORD69","SNORD71","SNORD74","SNORD76","SNORD8","SNORD81","SNORD82","SNORD84","SNORD86","SNORD87","SNORD89","SNORD94","SNORD96A","SNORD97","SNORD99","TMX1","ZFAS1","CDKN2B-AS1","CWF19L1","EIF4A1","EP400","HIF1A-AS2","MEG8","NOP56","RACK1","RPL13A","RPS13","SCARNA12","SNHG12","SNHG20","SNHG5","SNHG6","SNHG7","SNHG8","AP1G1","CFDP1","CHD8","DDX39B","EEF2","EIF4A2","EIF4G2","GAS5","GNL3","HSPA8","HSPA9","IPO7","MYRIP","NAN","NCL","NFATC3","PCAT4","PPAN","PRKAA1","PRRC2A","PTCD3","RABGGTB","RNF149","RPL10","RPL12","RPL13","RPL17","RPL21","RPL23","RPL23A","RPL4","RPL7A","RPLP2","RPS11","RPS2","RPS20","RPS3","RPS8","SF3B3","SLC25A3","SNHG1","SNHG3","SNHG9","SNORD1C","SNORD35B","SNORD37","SNRPB","SNX5","TAF1D","TNPO2","TOMM20","WDR43")
genes.in.dataset <- annotLookup[annotLookup$hgnc_symbol %in% gene.set.of.interest,]
genes.in.dataset.no.duplicates <- unique(genes.in.dataset$hgnc_symbol)


# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control","basal","her2","luminal a","luminal b"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE45827", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE45827", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("bottomleft", inset=c(0,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
#pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE45827")