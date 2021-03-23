GetAnnoFiles <- function(platform){
    platform <- as.character(platform)
    if(length(list.files(pattern=platform)) == 0){
        download.file(paste0("https://gemma.msl.ubc.ca/annots/", platform, "_noParents.an.txt.gz"), destfile=paste0(platform, ".gz"))
    }
    if(length(list.files(pattern=platform)) > 1){
        print("Multiple annotation files matching the platform exist")
    } else {
        warning("Using existing annotation file, consider updating")
    }
    Anno_file <- read.table(paste0(platform, ".gz"), comment="#", header=T, quote='"', sep="\t")
    return(Anno_file)
}

GetMarkers <- function(region = "Cortex"){
    region <- as.character(region)
    CellType_genes <- mouseMarkerGenesCombined[[region]]

    for (i in 1:length(CellType_genes)) {
        CellType_genes[[i]] <- as.vector(mouse2human(CellType_genes[[i]])$humanGene)
    }
    
    names(CellType_genes) <- sapply(names(CellType_genes), function(x) paste(x, "Genes", sep="_"))
    
    # Exclude specific genes: 5HTR3A - also expressed in VIP negative 5HTR3A cells.
    if (region == "Cortex") {
        CellType_genes$GabaVIPReln_Genes <- CellType_genes$GabaVIPReln_Genes[!CellType_genes$GabaVIPReln_Genes %in% "HTR3A"]
    }
    
    names(CellType_genes) <- sapply(names(CellType_genes), function(x) gsub(" ", "_", x))
    return(CellType_genes)
} 

import_and_filter <- function(files_per_cohort, col_names, tx2gene, gtf_path) {
 txi <- tximport::tximport(files_per_cohort, type = "salmon", tx2gene = tx2gene)
 cts <- txi$counts
 cts <- cts[rowSums(cts)>0,]
# rownames(cts) <- gsub("\\.\\d+", "", rownames(cts))
 txi$counts <- cts
 # add gene symbols to extract neuronal marker genes
 gene_symbols <- get_gene_info(g_source = "gff", gtf_path = gtf_path) %>% dplyr::select(gene_id, gene_name, gene_type, seqnames)
 counts_geneLevel <- as.data.frame(txi$counts)
 counts_geneLevel$gene_id <- rownames(counts_geneLevel)
# counts_geneLevel <- dplyr::left_join(counts_geneLevel, gene_symbols, by = "gene_id")

 ## FILTER 
 #filter mt genes
 # Get subset of mitochondrial and rRNA genes
 MitoGenes <- gene_symbols %>% dplyr::filter(seqnames == "chrM")
 RiboGenes <- gene_symbols %>% dplyr::filter(gene_type %in% c("rRNA", "Mt_rRNA"))
 ProtGenes <- gene_symbols %>% dplyr::filter(gene_type == "protein_coding")

 # Counts for mitochondrial genes
 MitoCountSum <- colSums(counts_geneLevel %>% 
			dplyr::filter(gene_id %in% MitoGenes$gene_id) %>%
			dplyr::select(-gene_id))
 # Counts excluding mitochondrial genes and mitochondrial reads fraction over
 # total lib. size
 MitoCountFiltered <- counts_geneLevel %>%
        	dplyr::filter(!(gene_id %in% MitoGenes$gene_id))
 MitoFiltCountSum <- colSums(MitoCountFiltered %>% dplyr::select(-gene_id))
 # For each sample, get the top 5 highest counted genes and normalise the counts
 # by the sample library size
 TopFiveProportionNoMT <- sapply(col_names, function(sid) { 
                  SubMatrix = MitoCountFiltered %>%
			  dplyr::select(gene_id, all_of(sid))
                  names(SubMatrix)[2] <- "Counts"
                  TopFive = SubMatrix %>%
			  dplyr::arrange(desc(Counts)) %>%
			  head(5)
                  TopFive %<>%  dplyr::mutate(Proportion = Counts/MitoFiltCountSum[sid])
                  Genes <- gene_symbols[match(TopFive$gene_id, gene_symbols$gene_id),] %>%
			  dplyr::select(-gene_id)
		  Genes$sample_id  <- sid
		  temp <- cbind(Genes, TopFive)
                  temp
		}, simplify = F) %>% data.table::rbindlist() %>% as.data.frame(stringsAsFactors=FALSE) 
	# Restrict to genes that gather at least 1% of the reads for the sample
	minProp <- 0.01
	TopFiveProportionNoMT %<>% dplyr::filter(Proportion > minProp) %>%
    		dplyr::mutate(gene_id = as.character(gene_id))
	# How much those top 5 genes represent for the sample in terms of coverage
	TopFiveSum <- TopFiveProportionNoMT %>%
    		dplyr::group_by(sample_id) %>%
    		dplyr::summarise(TotProp = sum(Proportion)) %>%
    		as.data.frame(stringsAsFactors=FALSE) %>%
    		dplyr::arrange(TotProp)
	TopFiveGeneFreq <- TopFiveProportionNoMT %>%
    		dplyr::group_by(gene_id) %>%
    		dplyr::summarise(n = n()) %>%
    		dplyr::ungroup(.) %>%
    		dplyr::mutate(prop = n / length(col_names)) %>%
    		dplyr::mutate(gene_id = as.character(gene_id),
			      ensemblID2 = paste0(gene_id, " (", n, "/", length(col_names), ")")) %>%
    		as.data.frame(stringsAsFactors=FALSE)
	TopFiveGeneFreq <- merge(TopFiveGeneFreq, gene_symbols[!duplicated(gene_symbols$gene_id),], by.x = "gene_id", by.y = "gene_id", all.x = T, all.y = F) %>%
    dplyr::arrange(desc(prop)) 

    	TopFiveTable <- TopFiveGeneFreq %>%
    		dplyr::select(gene_id, n, prop) %>%
    		droplevels %>%
    		dplyr::arrange(desc(prop))
  
	# Get the common genes with the highest count in all the samples (and the ribosomal gene..)
	CommonTopGenes <- TopFiveTable %>%
    		dplyr::filter(prop > 0.5)

	countMatrixFiltered <- MitoCountFiltered %>%
               dplyr::filter(!(gene_id %in% (CommonTopGenes$gene_id)))

	
 # filter low expressed
	GenesKept <- data.frame(matrix(rep(0, 3), ncol = 1))
	colnames(GenesKept) <- "Norway"
	rownames(GenesKept) <- c("Nuclear genes", "Expression above median", "With Gene Symbol")

	# Create log2 CPM matrix after removal of mitochondria-encoded genes
	cpmMatrixFiltered <- as.data.frame(log2(count2CPM(countMatrixFiltered %>% dplyr::select(-gene_id)) + 1))
	GenesKept[1,1] <- nrow(cpmMatrixFiltered)
	# Add gene symbols
	cpmMatrixFiltered %<>% dplyr::mutate(gene_id = countMatrixFiltered$gene_id)
	ExpDataCPM <- dplyr::left_join(cpmMatrixFiltered,
                                   gene_symbols	%>%
					   dplyr::select(gene_id, gene_name),
                                   by = "gene_id") %>%
                      dplyr::select(gene_id, gene_name, everything())

	# Calculate a cohort-specific noise level based on median CPM expression
	MaxNoise <- ExpDataCPM[, -c(1,2)] %>% as.matrix %>% median
	# ...
	MaxNoise  <- 0.1
	genesAboveNoise <- apply(ExpDataCPM[, -c(1,2)] %>% as.matrix, 1, function(Gene) {
                                quantile(Gene, 0.80) > MaxNoise
                              })

	# Filter genes below noise level
	ExpDataCPM <- ExpDataCPM[genesAboveNoise, ]
	GenesKept[2,] <- nrow(ExpDataCPM)
	# Filter genes with no Gene Symbol or with duplicated Gene Symbol
	ExpDataCPM <-   ExpDataCPM %>% 
				dplyr::filter(gene_name != "", !is.na(gene_name)) %>%
                             	dplyr::mutate(Sum = do.call(pmax, select_if(., is.numeric))) %>%
                             	dplyr::arrange(desc(Sum)) %>%
                             	dplyr::distinct(gene_name, .keep_all=TRUE) %>%
                             	dplyr::select(-Sum) %>%
                             	dplyr::arrange(gene_id)

	GenesKept[3,] <- nrow(ExpDataCPM)


	return(ExpDataCPM)
	}
get_est_gene_length <- function(files_per_cohort, col_names, tx2gene) {
 txi <- tximport::tximport(files_per_cohort, type = "salmon", tx2gene = tx2gene)
 gene_symbols <- get_gene_info(id_list = rownames(txi$length), tx = FALSE) %>% dplyr::select(gene_id,gene_name)
 print(colnames(txi$counts))
 tx_length <- as.data.frame(txi$length)
 colnames(tx_length)  <- col_names
 tx_length$gene_id <- rownames(tx_length)
 tx_length <- dplyr::left_join(tx_length, gene_symbols, by = "gene_id")
 return(tx_length)
}
get_len_scld_tpm <- function(files_per_cohort, col_names, tx2gene) {
 txi <- tximport::tximport(files_per_cohort, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM" )
 gene_symbols <- get_gene_info(id_list = rownames(txi$length), tx = FALSE) %>% dplyr::select(gene_id,gene_name)
# print(colnames(txi$counts))
 tx <- as.data.frame(txi$counts)
 colnames(tx)  <- col_names
 tx$gene_id <- rownames(tx)
 tx <- dplyr::left_join(tx, gene_symbols, by = "gene_id")
 return(tx)
}

pea <- function(res, gene_col = "Gene.names", sig_col = "P.Value", fc_col = "logFC", out_file = "./treemap_up_down.pdf", title = "") {
	res <- get_one_sided_pval(res, GeneCol = gene_col, pvalCol = sig_col, logFCcol = fc_col) 
        down <- ermineR::gsr(scores = res, scoreColumn="DownPval",
                                   bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                                   aspects=c("B", "M", "C"),
                                   iterations=200000)$results 
        up <- ermineR::gsr(scores  = res, scoreColumn="UpPval",
                                  bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                                   aspects=c("B", "M", "C"),
					 iterations=200000)$results 

# SIMPLIFYING PATHWAYS
library("colorspace")
maxP <- 200
    # UP
    obj <- up %>%
        filter(CorrectedPvalue < 0.05) %>%
        arrange(CorrectedPvalue) %>%
        select(Name, GeneMembers, CorrectedPvalue, NumGenes)
  #   print(head(obj))
    pathways <- apply(obj, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
    names(pathways) <- obj$Name
    pathways <- qdapTools::list2df(pathways, col1="gene", col2="Name")[,c(2,1)]
    clustersUP <- cluster_pathways(pathways, method="overlap", subsetsize=maxP)
    names(clustersUP) <- c("Name", "members")
# add fishers to aggr pvalues to sort clusters
    clustersUP <- left_join(clustersUP, obj %>% select(-GeneMembers, members=Name), by="members") %>%
        mutate(log10P=-log10(CorrectedPvalue)) %>%
        arrange(CorrectedPvalue)
    # DOWN
    obj <- down %>%
        filter(CorrectedPvalue < 0.05) %>%
        arrange(CorrectedPvalue) %>%
        select(Name, GeneMembers, CorrectedPvalue, NumGenes)
	pathways <- apply(obj, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
    names(pathways) <- obj$Name
    pathways <- qdapTools::list2df(pathways, col1="gene", col2="Name")[,c(2,1)]
    clustersDOWN <- cluster_pathways(pathways, method="overlap", subsetsize=maxP)
    names(clustersDOWN) <- c("Name", "members")
    clustersDOWN <- left_join(clustersDOWN, obj %>% select(-GeneMembers, members=Name), by="members") %>%
        mutate(log10P=-log10(CorrectedPvalue)) %>%
        arrange(CorrectedPvalue)


maxD <- 100
pdf(out_file)
    draw_treemap(dtf=bind_rows(clustersUP, clustersDOWN %>%
			       mutate(log10P=-log10P)) %>%
                         arrange(desc(abs(log10P))) %>%
                         head(n=maxD),
                 title=title)                               
dev.off()
return(list("clusters"= list("up" = up, "down" = down), "plot" = p))
}

rescale <- function(x) (x-min(x))/(max(x) - min(x)) 

transform_func <- function (m) {
	 log(m)
	 # m
	}

grepID <- function(string) {
		str_match(string,"[[:digit:]]+[[:punct:]][[:digit:]]+")
		    }

remove_zeros <- function(rawData, info, n = 6) {
#
# Now we can set to NA in the pData object all the samples from batches with
# non-detected proteins.
# Yes, this is a for loop :'-D
na_in_all_per_batch <- vector(mode = "list", length = length(unique(info$batch)))
names(na_in_all_per_batch) <- unique(info$batch)
for (batch in unique(info$batch)) {
    sample_ids <- info[info$batch==batch,]$reporter.intensity.id
    samples_in_batch <- rawData %>% dplyr::select(all_of(sample_ids))
    row_na <- which(rowSums(samples_in_batch)==0) 
    na_in_all_per_batch[[batch]] <- row_na    
}
#aThird<-(length(info$sample_id)*(2/3))
lapply(na_in_all_per_batch,length)
# which proteins are zero in either more than 7 batches
no_of_batches <- table(unlist(na_in_all_per_batch))
allZero_in_over_n_batches <-names(no_of_batches)[which(no_of_batches > n)]

## Find rows with no variance
noVarProts <-  which(apply(rawData %>% dplyr::select( -Gene.names, -Protein.IDs), 1, var)==0)
removedZero <- rawData[-c(unique(noVarProts, allZero_in_over_n_batches)),]
}
loading_norm <- function(df, info, within_batch = TRUE, use_median = TRUE, name_cols = c("Gene.names","Protein.IDs"), method = "lib_size") {
if (within_batch == TRUE) {
 orig_order <- colnames(df %>% dplyr::select(-name_cols))
 normDF <- df %>% dplyr::select(name_cols)
# if(use_median == TRUE) {
#  avg_lib_size <- median(colSums(df %>% dplyr::select(-name_cols)))
# } else {
#   avg_lib_size <- mean(colSums(df %>% dplyr::select(-name_cols)))
# }
 for (batch in unique(info$batch)) {
    batchCols <- info[info$batch == batch, "reporter.intensity.id"]
    m <- as.matrix(dplyr::select(df, all_of(batchCols)))
    # norm_vec <- log(avg_lib_size) / log(colSums(m))
    # m_norm <- sweep(m, 2, norm_vec, FUN = "*")
    if (method == "log_lib_size") {
     m_norm <- apply(m, 2, function(x) { log(x) / log(sum(x))})
    } else if (method == "lib_size") {
     m_norm <- apply(m, 2, function(x) { x / sum(x)})
    }  else if (method == "quantile") {
     m_norm <- normalize.quantiles.robust(m, use.median = use_median, remove.extreme = "both") 	
    } else {
	    stop("Available methods: lib_size, log_lib_size, quantile")
    }
    batch_df_norm <- as.data.frame(m_norm)
    cnames <- colnames(normDF)
    normDF <- cbind(normDF, batch_df_norm)
    colnames(normDF)<-c(cnames, batchCols)
 }
 normDF %<>% dplyr::select(all_of(colnames(df)))
} else {
 normMat <-  normalize.quantiles.robust(as.matrix(df[, -c(1, 2)]), remove.extreme = "both", use.median = use_median)
 normDF <- cbind(df[,c(1,2)],normMat)
 colnames(normDF)<- colnames(df)
# normDF  <- NULL
# print("need to code this")
 }
 return(normDF)
}

batch_corr <- function(rawData, info, name_cols = c("Gene.names", "Protein.IDs")) {
 # reporter.intensity.id of all the references of each channel (channel 0, batch 1:10)
 ref_ids <- info %>% dplyr::filter(condition1 == "Reference") 
 referenceRawCounts <- rawData %>% dplyr::select(ref_ids$reporter.intensity.id, all_of(name_cols))
 # calculate correction factor
 cf <- referenceRawCounts %>% dplyr::mutate(denom = rowSums(.[1:10]) / 10) %>%
	 mutate_at(grep("Reporter.intensity", colnames(.)), ~ . / (denom))
 # rename refchannel cf to batch index
 cfBatchIdx <-colnames(cf)[colnames(cf) %>% grep("Reporter.intensity", .)] %>%
	 stringr::str_match(.,"[:digit:]+$")
 colnames(cf)[colnames(cf) %>% grep("Reporter.intensity",.)]<-cfBatchIdx
	
 # apply a function to all measurement columns
 #extract the column name of the column given by dplyr
 # extract only the batch index from it with string match (the last digit (one or more digit)
 #if either the rawData value or the cf value (cf value should be from the corresponding protein and batch) is zero 
 #dividing by cf will result in either nan or infinity, in this case do not divide by cf --> THINK ABOUT THAT
 rawData %>% dplyr::mutate_at(grep("Reporter.intensity",colnames(.)), funs({
			colName <- stringr::str_match(quo_name(quo(.)),"[:digit:]+$")
			cf_val  <- cf[,colName]
			 # print(.) 
			 # print(paste0("intensity: ",.)) 
			 # print(paste0("cf_val: ",cf_val))
		 	 ifelse((is.infinite(./cf_val) | is.nan(./cf_val)),.,(./cf_val))				
		#	. / (1 + cf_val)
	})) -> normCN 
# print(apply(cf,2,function(x)(sum(is.infinite(x)))))
 #p <- ggplot(reshape2::melt(cf), aes(x = value)) + 
#	geom_histogram() +
#	facet_wrap(~ variable, scales = "free")
# ggsave(plot = p, filename = "./cf_dist.pdf", device = "pdf")
 return(normCN)
}

proteus_norm <- function(evidenceFile,info,rawData, normalize = T, protein_col = "protein_group") {
 measCols <-  c(paste0("Reporter intensity ",unique(info$channel)))
 names(measCols) <- gsub(" ",".",measCols)
 evi <- proteus::readEvidenceFile(evidenceFile, measure.cols = measCols, data.cols = c(evidenceColumns,"Gene names"))
 colnames(evi)[9] <- "Gene.names"
 metaProteus <- info %>% 
  dplyr::mutate(condition2 = ifelse(age_years < 1 | age_years > 60 | is.na(age_years), condition1, "MD_Control")) %>%  
  dplyr::mutate(condition2 = ifelse(condition2 == "Control" & cohort == "Netherlands Brain Brank", Control_NBB, condition2)) %>%
  dplyr::rename(experiment = batch, old_condition = condition, condition = condition2 ) %>%
  dplyr::mutate(sample=paste0("Reporter.intensity.",channel,".",experiment),
		measure=paste0("Reporter intensity ",channel))  
 pepdat <- makePeptideTable(evi, metaProteus, protein.col = protein_col, measure.cols = measCols, aggregate.fun = aggregateMedian, experiment.type = "TMT")
 # create protein table from evidence file and aggregate peptides
 prodat <- makeProteinTable(pepdat, aggregate.fun = aggregateHifly, hifly = 3)
 m <- prodat$tab %>% as.data.frame(.) %>%
	 dplyr::mutate(Protein.IDs = rownames(.))
 rownames(m) <- rownames(prodat$tab)
 if(normalize == TRUE) {
  # normalize WITHIN batch with CONSTANd algorithm
  prodat <- normalizeTMT(prodat)
  m <- prodat$tab %>% as.data.frame(.) %>%
	  dplyr::mutate(Protein.IDs = rownames(.))
  rownames(m) <- m$Protein.IDs
 }
 if (protein_col == "protein_group") {
  m <- m[rawData$Protein.IDs,] %>% 
	dplyr::left_join(., rawData[,c("Gene.names", "Protein.IDs")], by = "Protein.IDs") %>%
	dplyr::select(all_of(colnames(rawData)))
 } else {
  m <- as.data.frame(prodat$tab) %>% dplyr::mutate(Protein.IDs = rownames(.))
  annot <- fetchFromUniProt(prodat$proteins, columns = c("database(Ensembl)", "genes(PREFERRED)"))
  colnames(annot) <- c("Protein.IDs", "tx_id", "Gene.names")
  # annot %>% separate(tx_id, c("alt_tx_id", "tx_id"), sep = ";") %>%
  # dplyr::left_join(., DTU::get_gene_info(tx_id, tx = T)[,c("gene_id","tx_id")], by = "tx_id")
  m <- dplyr::left_join(m, annot, by = "Protein.IDs")
 }
 m %<>% dplyr::mutate_if(is.numeric, funs(tidyr::replace_na(., 0)))
 return(list(m = m, obj = prodat))
  
}
tryQuery <- function(url, maxtry=5) {
  for(i in 1:maxtry) {
    res <- tryCatch(read.delim(url, stringsAsFactors=FALSE), error=function(err) NULL)
    if(!is.null(res)) return(res)
    Sys.sleep(3)
  }
  stop(paste("UniProt not responding."))
}
fetchFromUniProt  <- function (unis, columns = c("genes", "protein names"), col.names = NULL,
    batchsize = 400, verbose = FALSE)
{
    # good.cols <- columns %in% allowedUniProtColumns
    #if (sum(!good.cols) > 0) {
     #   bad.cols <- paste(columns[!good.cols], collapse = ",")
      #  stop(paste("Columns", bad.cols, "are not allowed in UniProt queries. Here is a list of allowed columns:",
       #     paste(allowedUniProtColumns, collapse = ", "), "."))
    #}
    if (is.null(col.names)) {
        col.names <- columns
    }
    else {
        if (length(col.names) != length(columns))
            stop("col.names has to be the same length as columns.")
    }
    columns <- c("id", columns)
    col.names <- c("id", col.names)
    colstr <- paste(columns, collapse = ",")
    unis <- as.character(na.omit(unis))
    acc.test <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",
        unis, perl = TRUE)
    unis <- unis[acc.test]
    if (sum(!acc.test) > 0) {
        warning("Some identifiers do not conform to UniProt accession number format. Skipping.")
    }
    if (length(unis) == 0) {
        stop("No valid UniProt accession numbers found.")
    }
    url <- "http://www.uniprot.org/uniprot/"
    batches <- split(unis, ceiling(seq_along(unis)/batchsize))
    nbatch <- length(batches)
    res <- lapply(seq_along(batches), function(i) {
        if (verbose)
            cat(paste("Batch", i, "out of", nbatch, "\n"))
        ids <- batches[[i]]
        qry <- paste(paste0("id:", ids), collapse = "+or+")
        qurl <- URLencode(paste0(url, "?query=", qry, "&format=tab&columns=",
            colstr))
        tryQuery(qurl)
    })
    df <- do.call(rbind, res)
    colnames(df) <- col.names
    df
}

gsea_stringdb <- function(obj, res) {
 library(STRINGdb)
 string_db <- STRINGdb$new(species = 9606, score_threshold = 0, version = "10")
 # define gene background 
 bg <- unique(gsub(";.*","",obj$proteins))
 bg <- string_db$map( data.frame(protein=bg), "protein", removeUnmappedRows = T )   
 # Load stringdb data for homo sapiens
 string_db <- STRINGdb$new(species = 9606, score_threshold = 0, version = "10", backgroundV = bg$STRING_id)
 # define hits
 hits <- unique(gsub(";.*","",subset(res,significant == TRUE)$protein))
 #map gene names to what is needed for stringdb
 hits <- string_db$map( data.frame(protein=hits), "protein", removeUnmappedRows = T )   
 # "category"      category for which to compute the enrichment (i.e. "Process", "Component", "Function", "KEGG", "Pfam", "InterPro", "Tissue", "Disease".  default "Process")
 bp <- string_db$get_enrichment(hits$STRING_id, category = "Process",
	methodMT = "fdr", iea = TRUE ) %>%
	dplyr::filter(pvalue_fdr < 0.05) %>%
	tbl_df #%>% print(n = Inf)
 cc <- string_db$get_enrichment(hits$STRING_id, category = "Component",
	methodMT = "fdr", iea = TRUE ) %>%
	dplyr::filter(pvalue_fdr < 0.05) %>%
	tbl_df #%>% print(n = Inf)
 mf <- string_db$get_enrichment(hits$STRING_id,
	category = "Function", methodMT = "fdr", iea = TRUE ) %>%
	dplyr::filter(pvalue_fdr < 0.05) %>%
	tbl_df #%>% print(n = Inf)
return(list("bp" = bp, "cc" = cc, "mf" = mf, "stringdb" = string_db))
}

# enriched <- gsea_stringdb(obj,res)
# little helper function to create the input for clustering of pathways
create_input <- function(enriched, res, categoryTmp) {
 stringDBRes <- enriched[[categoryTmp]]
 pathways <- stringDBRes$pvalue_fdr
 names(pathways) <- stringDBRes$term_description
 res <- na.omit(res %>% dplyr::select(adj.P.Val, protein) %>%
	dplyr::mutate(protein = gsub(";.*","",protein)))
 # load all annotations and subset by category of intrest e.g. "Process"
 ann = enriched[["stringdb"]]$get_annotations() %>%
	dplyr::filter(category == dplyr::case_when(
		categoryTmp == "bp" ~ "Process",
		categoryTmp == "cc" ~ "Component",
		categoryTmp == "mf" ~ "Function")) %>%
        dplyr::filter(term_id %in% stringDBRes$term_id)
 # create mapping to gene name
 id_map <- enriched[["stringdb"]]$map(data.frame(protein = res$protein), "protein", removeUnmappedRows = T)   
 ann %<>% dplyr::left_join(., id_map, by = "STRING_id")
 ann <- na.omit(ann) # basically i dont have all the genes in my data that are in the annotation so i reduce it to my background
 ann %<>% dplyr::left_join(., res ,by = "protein")
 info <- stringDBRes %>% dplyr::select(term_id, term_description)
 ann %<>% left_join(., info, by = "term_id") #%>% dplyr::select(protein, STRING_id, adj.P.Val) 
 psets <- split(ann,f = ann$term_description)
 psets <- psets[names(pathways)]
 psets <- lapply(psets, function(set) {
	set %<>% dplyr::select(protein, STRING_id, adj.P.Val)
	colnames(set) <- c("gene", "gene_id", "pvalue")
        return(set)
 })
 p_obj <- list(pathways = as.list(pathways), genesets = psets)
 return(p_obj)
 }
count2CPM <- function(countData){
  apply(countData, 2, function(SampleCount){
    SampleLibSize = sum(SampleCount)
    (10^6)*(SampleCount/SampleLibSize)
  })
}
investigateTopGenes <-function(countMat, name_cols = c("gene_id", "gene_name"))
{
	top <- apply(countMat %>% dplyr::select(-c(name_cols)),2,function(col){
					      	top = which.max(col)
						for (i in 1:4)
						{	col[top] = 0
							top <- c(top,which.max(col))
						}	
	    								return(top)
	     })
#	print(head(top))
gene_idcs <- as.numeric(as.vector(names(sort(table(top), decreasing = T))))
topGenesInAllSamples <- data.frame(countMat[gene_idcs, name_cols],
				   gene_idx = gene_idcs,
				   top5InNoSampl = as.vector(sort(unname(table(top)), decreasing = T)))
countMat %<>% dplyr::select(-name_cols)
propReadsInTopGenesPerSampl <- sapply(seq(1,ncol(top)),function(i) {
#	         print(countMat[top[, i],i])
#		 print(sum(countMat[,i]))
	         perc = sum(countMat[top[, i],i])/sum(countMat[,i]) * 100
		return(perc)
	     })
	     
return(list(topGenesInAllSamples, propReadsInTopGenesPerSampl,top))
}
k_plot <- function(mydata,kMax) {
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:kMax) wss[i] <- sum(kmeans(mydata,nstart=25,iter.max=10000,centers=i)$withinss)
plot(1:kMax, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
}
expVar <- function(pca) { pca$sdev^2 / sum(pca$sdev^2)}
rescale <- function(x) (x-min(x))/(max(x) - min(x)) 


DESeq2RUN <- function(data, Meta, model){
    DESeqDS <- DESeqDataSetFromTximport(txi = data, colData = Meta, design = model)
    DESeqOut <- DESeq(DESeqDS)
    return(DESeqOut)
}
GetDESeq2Results <- function(DESeqOut, coef, alpha=0.05, indepFilter=TRUE, geneNames=NULL){
#    if (is.null(geneNames)) {
#        geneNames <- biomaRt::getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "gene_biotype", "chromosome_name"), 
#                           mart=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
#    }
#    DEresults$GeneSymbol <- geneNames$hgnc_symbol[match(rownames(DEresults), geneNames$ensembl_gene_id)]
    DEresults$GeneSymbol <- rownames(DEresults)
#    DEresults %<>% data.frame #%>% dplyr::filter(GeneSymbol != "")
    return(DEresults)
}
## perform pca on randomly sampled proteins
run_pca <- function(data, indices) {
 d<- data[indices,] # allows boot to select sample
 pca <- prcomp(t(d),scale=T)
 pc1 <- pca$x[,1]
# corVal <- cor(pc1,as.numeric(as.vector(info$deficiancy)))
 expVar <- pca$sdev^2
 expVar <- expVar[1] / sum(expVar)
# euklInfo <- distCentrPCA(pca,info$condition)
# cos <- pca
# return(c("expVar"=expVar,"cor"=corVal,"distGroups"=euklInfo[[1]],"centroids"=list(euklInfo[[1]],euklInfo[[2]])))
return(c(exp_var=expVar))
}
cluster_all_GO <- function(obj, res) {
 categories <- c("mf", "cc", "bp")
 enriched_obj <- gsea_stringdb(obj, res) 
 p_obj_lst<- lapply(categories, function(cn){
  create_input(enriched_obj, res, cn)
 })
 names(p_obj_lst) <- categories

 # cluster pathways
 clusters <- lapply(seq(1,length(p_obj_lst)), function(i) {
  run_pathway_vis(input = p_obj_lst[[i]], plotTitle = categories[[i]],
	 subsetsize = 50, outPDF = paste0(plotOutDir, names(p_obj_lst)[i],".pdf"))
 })
 names(clusters) <- categories
 return(clusters)
}
draw_treemap <- function(dtf, title="", vp=NULL) {
    treemap::treemap(
        dtf=dtf,
        index=c("Name", "members"), 
        vSize="NumGenes", 
        vColor="log10P",
        sortID="color",
        title=title,
        type="value",
        palette=rev(brewer.pal(10, "RdBu")),
        border.col=c("white", "black"),
        border.lwds=c(5,1),
        overlap.labels=1,
        bg.labels="#FFFFFFCC",
        force.print.labels=FALSE,
        vp=vp)
}

gene_pathway_matrix <- function(input){
    m <- table(input[[2]], input[[1]]) == 1
    m[,match(unique(input[[1]]), colnames(m))]
}
sim_matrix <- function (gmmat){
    #require(parallel)
    Kappas <- sapply(1:ncol(gmmat), function(p1) {
                  sapply(1:ncol(gmmat), function(p2) irr::kappa2(gmmat[,c(p1,p2)])$value)
    })
    rownames(Kappas) <- colnames(gmmat)
    colnames(Kappas) <- colnames(gmmat)
    return(Kappas)
}
jaccard_matrix <- function (gmmat){
    as.matrix(1 - dist(t(gmmat), method="binary"))
}
overlap_matrix <- function (gmmat){
    #require(parallel)
    Overlap <- sapply(1:ncol(gmmat), function(p1) {
                   sapply(1:ncol(gmmat), function(p2) {
                       sum(gmmat[,p1] & gmmat[,p2]) / min(sum(gmmat[,p1]), sum(gmmat[,p2]))
                   })
               })
    rownames(Overlap) <- colnames(gmmat)
    colnames(Overlap) <- colnames(gmmat)
    return(Overlap)
}
update_simmat <- function(simmat, gmmat, failed, passed){
    diag(simmat) <- 1
    #gene_pathway_mat <- gene_pathway_matrix(pathways)
    genes_in_cluster <- rep(0, length(rownames(gmmat)))
    names(genes_in_cluster) <- rownames(gmmat)
    # genes in the union of the two pathways
    merged <- unique(c(rownames(gmmat)[gmmat[,failed]],
                       rownames(gmmat)[gmmat[,passed]]))
    genes_in_cluster <- ifelse(names(genes_in_cluster) %in% merged, 1, 0)
    names(genes_in_cluster) <- rownames(gmmat)
    new_kappa_col <- sapply(1:nrow(simmat), function(x){
            irr::kappa2(data.frame(gmmat[,x],genes_in_cluster))$value
    })
    simmat[,passed] <- new_kappa_col
    simmat <- simmat[-which(rownames(simmat)==failed), -which(rownames(simmat)==failed)]
    diag(simmat) <- NA
    return(simmat)
}
cluster_pathways <- function(pathways, threshold=0.4, subsetsize=length(unique(pathways[[1]])),
                            lowThreshold=50, highThreshold=1000, method=c("kappa", "jaccard", "overlap")) {
    # check pathways input is data.frame
    stopifnot("data.frame" %in% class(pathways))
    # max number of pathways allowed are 1000, safety measure
    totPaths <- length(unique(pathways[[1]]))
    if (subsetsize > totPaths) {
        warning(paste0("Subset size chosen (", subsetsize, ") greater than the number of pathways (", totPaths, "), setting subsetsize <- ", totPaths))
    	subsetsize <- totPaths
    }
    if (subsetsize > 1000){
        warning("Maximum subset size allowed is 1000, setting subsetsize <- 1000")
    	subsetsize <- 1000
    }

    method <- method[1]
    
    # Filter to only top pathways
    pathways[,1] <- as.character(pathways[,1])
    pathways[,2] <- as.character(pathways[,2])
    pathways <- pathways[pathways[[1]] %in% head(unique(pathways[[1]]), n=subsetsize),]

    # membership of genes to pathways and pathway sizes
    membership_matrix <- gene_pathway_matrix(pathways)
    pSizes <- colSums(membership_matrix)

    # similarity of pathways
    if (method == "kappa") {
        simmat <- sim_matrix(membership_matrix)
    } else if (method == "jaccard") {
        simmat <- jaccard_matrix(membership_matrix)
    } else if (method == "overlap") {
        simmat <- overlap_matrix(membership_matrix)
    } else {
        stop(paste0("Method \"", method, "\" unknown"))
    }
    diag(simmat) <- NA

    pathwayNames <- colnames(simmat)   # pathway names
    if (any(pathwayNames != unique(pathways[[1]]))) stop("HORRIBLE ERROR")

    cluster <- matrix(rep(0, subsetsize*subsetsize), nrow=subsetsize)
    diag(cluster) <- 1
    colnames(cluster) <- pathwayNames
    rownames(cluster) <- pathwayNames

    iter <- 1
    while(max(simmat[upper.tri(simmat, diag=FALSE)]) >= threshold){
        message(paste0("Iteration ", iter, "..."))
        index <- which(simmat == max(simmat[upper.tri(simmat, diag = FALSE)]), arr.ind=TRUE)
        t_i = index[1,1] # first row index
        t_j = index[1,2] # first col index
        t_both <- c(t_i, t_j)
        pair_size <- c(pSizes[rownames(simmat)[t_i]], pSizes[rownames(simmat)[t_j]])
        pair_names <- c(rownames(simmat)[t_i], rownames(simmat)[t_j])
        if ( (pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold) && (pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold) ) {
            failed <- rownames(simmat)[max(t_both)]
            passed <- rownames(simmat)[min(t_both)]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold ) {
            failed <- rownames(simmat)[t_i]
            passed <- rownames(simmat)[t_j]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold ) {
            failed <- rownames(simmat)[t_j]
            passed <- rownames(simmat)[t_i]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else {
            failed <- rownames(simmat)[max(t_both)]
            passed <- rownames(simmat)[min(t_both)]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        }
        iter <- iter + 1
    }

    # remodel "cluster" structure to extract the titles
    cluster_idx <- which(diag(cluster)==1)
    clusters <- apply(cluster[,cluster_idx], 2, function(x) names(which(x==1)))
    data.frame(title=rep(names(clusters), sapply(clusters, length)),
               member=unlist(clusters, , FALSE),
               stringsAsFactors=FALSE,
               row.names=NULL)
}

get_one_sided_pval <- function(ResultsObj, adjust="BH", logFCcol="logFC", GeneCol="Gene.names", pvalCol="P.Value"){
    DESeqResultsDF <- data.frame(ResultsObj)
    #Just for now - remove the duplicated genes (5 at this point) 
    DESeqResultsDF <- DESeqResultsDF[!duplicated(DESeqResultsDF[,GeneCol]),]
    DESeqResultsDF$DownPval <- apply(DESeqResultsDF %>% dplyr::select(logFCcol, pvalCol), 1, function(x){
        if(x[1] < 0){
            x[2]/2
        } else {
            1-x[2]/2
        }
    })
    DESeqResultsDF$DownPvalAdj <- p.adjust(DESeqResultsDF$DownPval, "BH")
    DESeqResultsDF$UpPval <- apply(DESeqResultsDF %>% dplyr::select(logFCcol, pvalCol), 1, function(x){
        if(x[1] > 0){
            x[2]/2
        } else {
            1-x[2]/2
        }
    })
    DESeqResultsDF$UpPvalAdj <- p.adjust(DESeqResultsDF$UpPval, "BH")
    rownames(DESeqResultsDF) <- as.character(DESeqResultsDF[[GeneCol]])
    return(DESeqResultsDF)
}

GetAnnoFiles <- function(platform){
    platform <- as.character(platform)
    if(length(list.files(pattern=platform)) == 0){
        download.file(paste0("https://gemma.msl.ubc.ca/annots/", platform, "_noParents.an.txt.gz"), destfile=paste0(platform, ".gz"))
    }
    if(length(list.files(pattern=platform)) > 1){
        print("Multiple annotation files matching the platform exist")
    } else {
        warning("Using existing annotation file, consider updating")
    }
    Anno_file <- read.table(paste0(platform, ".gz"), comment="#", header=T, quote='"', sep="\t")
    return(Anno_file)
}

#' create annotation df
#' 
#' Create transcript annotation file to map transcripts to genes. Either with ensembldb or read from supplied gtf path. Used gtf file in referenceData dir.
#' @param version string ensembldb version (just the number). Default=75
#' @param g_scource character string genome annotation source, ensembldb or ucsc (->gtf file). Remember there were some modifications needed to be done for the gtf file.
#' @param ntx_filter integer threshold for maximum number of transcripts per gene. Genes with more transcripts will be removed from the annotation 
#' @return annotation dataframe as required by DRIMSeq in DTU::prep_dtu() to build DRIMSeqData object
#' @importFrom rlang .data 
#' @export
get_annot <- function(version = "86", g_source = "ucsc", gtf_path = NULL, ntx_filter = NULL, strip_version = T) {
 if (g_source == "ensembldb") {
  db <- paste0("EnsDb.Hsapiens.v", version)
  if (!(db %in% rownames(utils::installed.packages()))) {
   BiocManager::install(db)
  }
  if (!("ensembldb" %in% rownames(utils::installed.packages()))) {
   BiocManager::install("ensembldb")
  }
  library(db, character.only = TRUE)
  # Ensembl mapping between genes and transcripts
  edb <- get(db)
  txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene",
   columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid")))
  txdb <- subset(txdb, ensembldb::seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
  txdb <- txdb[, c("tx_id", "gene_id", "gene_name", "tx_biotype", "entrezid")]
  names(txdb) <- c("TXNAME", "GENEID", "GENENAME", "BIOTYPE", "ENTREZ")
  annot <- txdb[, c(2, 1)]
  tab <- table(annot$GENEID)
  annot$ntx <- tab[match(annot$GENEID, names(tab))]
  if (!(is.null(ntx_filter))) {
   annot <- subset(annot, ntx <= ntx_filter)
  }
  cat(crayon::blue("ensembl db obj info"))
  print(edb)
  return(annot)
 }
 else if (g_source == "ucsc") {
  if (is.null(gtf_path)) {
   cat(crayon::blue("you chose to load annotation from ucsc gtf file, you didnt provide a path to your gtf file (e.g. ucsc.hg19.gtf)"))
   return(NULL)
  } else {
   txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
   #annot <- ensembldb::select(txdb, ensembldb::keys(txdb, "GENEID"), "TXNAME", "GENEID", "GENENAME")

	# Create the tx2gene
	annot <- mapIds(txdb,
                  keys=keys(txdb, "GENEID"),
                  column="TXNAME",
                  keytype="GENEID",
                  multiVals="list") %>%
  	enframe(name="GENEID", value="TXNAME") %>%
  	unnest(cols=c("TXNAME")) %>%
  	dplyr::select("GENEID", "TXNAME") %>%
  	distinct()

      	tab <- table(annot$GENEID)
   	annot$ntx <- tab[match(annot$GENEID, names(tab))]
   if (!(is.null(ntx_filter))) {
    annot <- subset(annot, ntx <= ntx_filter)
    }
   if (strip_version == TRUE) {
   annot %<>% dplyr::mutate(GENEID = gsub("\\.\\d+", "", GENEID),
			     TXNAME = gsub("\\.\\d+", "", TXNAME))
   }
   return(annot)
  }
 } else {
  cat(crayon::blue("The database you want to use does not exist or not offered to be used within this function"))
  cat(crayon::blue("Will return NULL"))
  return(NULL)
 }
}



#' Retrieve gene annotation from ensembldb
#' 
#' Especially used for making the paper figures. To retrieve transcript biotypes and gene names based on either transcript id or gene id.
#' @param id_list character vector of ids used to retrieve gene annotation, if tx==T then it must be ensembl transcript ids that are valid in the ensembl release that is specified with the parameter version. Else it must be ensembl gene ids.
#' @param version string specifying ensembl release version (default is 75 which was used within this analysis, the last of grch37)
#' @param g_source char string specifing the database used to retrieve annotation (this is now redundant, removed functionality to use ucsc)
#' @param tx bool specifying whether the vector provided with id_list are transcript (tx==TRUE) or gene ids
#' @param exons
#' @param uniprot
#' @param symbol
#' @param promoters
#' @param upBP
#' @param downBPi
#' @return dataframe of requested gene annotation.
#' @importFrom dplyr %>%
#' @export
get_gene_info <- function(id_list, version = "86", g_source = "gff", gtf_path = NA, tx = FALSE, exons=F, uniprot = FALSE, symbol = FALSE, promoters = FALSE, upBP=2000, downBP=400) {
  if(g_source == "gff") {
	stopifnot(file.exists(gtf_path)) 
	# Create mapping to gene types
	gene_names <- rtracklayer::import(gtf_path) %>% as_tibble %>%
    	dplyr::select(seqnames, gene_id, gene_name, transcript_id, gene_type) %>%
    	na.omit %>%
   	distinct(gene_id, .keep_all=TRUE)
     return(gene_names)

  }
  if (g_source == "ensembldb") {
  db <- paste0("EnsDb.Hsapiens.v", version)   
  if (!(db %in% rownames(utils::installed.packages()))) {
   BiocManager::install(db)
  }
  if (!("ensembldb" %in% rownames(utils::installed.packages()))) {
   BiocManager::install("ensembldb") 
  }
  #not sure how to solve the porblem of having a library call in here
  library(db, character.only = TRUE)
  # Ensembl mapping between genes and transcripts
  edb <- get(db)
  if (symbol == TRUE) {
   print("Using GeneNameFilter")
   genes <- as.data.frame(ensembldb::genes(edb, filter = AnnotationFilter::GeneNameFilter(id_list), columns = c("gene_id", "gene_name", "gene_biotype", "entrezid")))
   return(genes)
  }
  if (tx == FALSE) {
   genes_gr <- ensembldb::genes(edb, filter = AnnotationFilter::GeneIdFilter(id_list), columns = c("gene_id", "gene_name", "gene_biotype", "entrezid"))
   genes <- as.data.frame(genes_gr)
   genes <- subset(genes, seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
   if (promoters == TRUE) {
    return(promoters(genes_gr, upstream = upBP, downstream = downBP))	
   } else {
    return(genes)
   }
  } else {
   if (exons == FALSE) {
    if (uniprot == TRUE) {
     txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid", "tx_name", "protein_id", "uniprot_id")))
     txdb <- txdb %>% dplyr::group_by(.data$tx_id) %>% dplyr::mutate(uniprotID = as.character(paste(.data$uniprot_id, collapse = ","))) %>% dplyr::select(-c(.data$uniprot_id, .data$entrezid)) %>% dplyr::distinct()
    } else {
     txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid", "tx_name")))
    }
    txdb <- subset(txdb, seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
    tab <- table(txdb$gene_id)
    txdb$ntx <- tab[match(txdb$gene_id, names(tab))]
   } else {
    if (!(tx == TRUE)) {
     print("please enter tx_ids and put tx==T to get exons as they are retrieved from ensembldb by tx_id")
    } else {
     txdb <- as.data.frame(ensembldb::exonsBy(edb, by = "tx", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "exon_id", "exon_idx", "exon_seq_start", "exon_seq_end", "seq_strand")))
     txdb <- subset(txdb, seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
     txdb <- txdb %>% dplyr::group_by(.data$tx_id)
    }    
   }
   return(txdb)
  }
 }
}
