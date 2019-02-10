#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tximport))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(rjson))
suppressPackageStartupMessages(require(workflowscriptscommon))

usage <- "tximport.R --files kallisto_results.txt --tx2gene tx2gene --type software_type --output-file output_file  [options]"

option_list <- list(
  make_option(c("-f", "--files"), type='character', dest='files', default=NA, help="Text file containing a list of filenames for the transcript-level abundances (one per line). Sample names will be derived from directory names"),
  make_option(c("-t", "--type"), type='character', dest='type', default='none', help="The type of software used to generate the abundances. Must be one of 'none', 'salmon', 'sailfish', 'kallisto', 'rsem', 'stringtie'. This argument is used to autofill the arguments below (geneIdCol, etc.) 'none' means that the user will specify these columns.", action='callback', callback=wsc_choose_from, callback_args=list(choices=c('none', 'salmon', 'sailfish', 'kallisto', 'rsem', 'stringtie'))),
  make_option(c("-o", "--outputCountsFile"), type='character', dest='outputCountsFile', default=NA, help="Counts output file path. Where output format is 'sparse', this should be a directory path"),
  make_option(c("-q", "--outputAbundancesFile"), type='character', dest='outputAbundancesFile', default=NA, help="Abundances output file path. Where output format is 'sparse', this should be a directory path"),
  make_option(c("-w", "--outputStatsFile"), type='character', dest='outputStatsFile', default=NA, help="File in which to output a summary of statistics (Kallisto only)"),
  make_option(c("-s", "--outputFormat"), type='character', dest='outputFormat', default='sparse', help="Output file format. Once of 'tsv' (tab separated), 'sparse' (Cellranger Matrix Market format), 'HDF5' (HDF5)", action='callback', callback=wsc_choose_from, callback_args=list(choices=c('tsv', 'sparse', 'HDF5'))),
  make_option(c("-x", "--txIn"), type='logical', dest='txIn', default=TRUE, help="Whether the incoming files are transcript level (default TRUE)"),
  make_option(c("-y", "--txOut"), type='logical', dest='txOut', default=FALSE, help="Whether the function should just output transcript-level (default FALSE)"),
  make_option(c("-c", "--countsFromAbundance"), type='character', dest='countsFromAbundance', default='no', help="Whether to generate estimated counts using abundance estimates: 'no' (default), 'scaledTPM' (scaled up to library size), 'lengthScaledTPM' (scaled using the average transcript length over samples and then the library size), or 'dtuScaledTPM' (scaled using the median transcript length among isoforms of a gene, and then the library size (dtuScaledTPM))", action='callback', callback=wsc_choose_from, callback_args=list(choices=c('no', 'scaledTPM', 'lengthScaledTPM', 'dtuScaledTPM'))),
  make_option(c("-g", "--tx2gene"), type='character', dest='tx2gene', default=NULL, help="A two-column tab-delimited text file linking transcript id (column 1) to gene id (column 2). This argument is required for gene-level summarization for methods that provides transcript-level estimates only (kallisto, Salmon, Sailfish)"),
  make_option(c("-v", "--varReduce"), type='logical', dest='varReduce', default=FALSE, help="whether to reduce per-sample inferential replicates information into a matrix of sample variances variance (default FALSE)"),
  make_option(c("-d", "--dropInfReps"), type='logical', dest='dropInfReps', default=FALSE, help="Whether to skip reading in inferential replicates (default FALSE)"),
  make_option(c("-i", "--ignoreTxVersion"), type='logical', dest='ignoreTxVersion', default=FALSE, help="logical, whether to split the tx id on the '.' character to remove version information, for easier matching with the tx id in gene2tx (default FALSE)"),
  make_option(c("-n", "--ignoreAfterBar"), type='logical', dest='ignoreAfterBar', default=FALSE, help="logical, whether to split the tx id on the '|' character (default FALSE)"),
  make_option(c("-l", "--geneIdCol"), type='character', dest='geneIdCol', default=NA, help="Name of column with gene id. if missing, the gene2tx argument can be used"),
  make_option(c("-m", "--txIdCol"), type='character', dest='txIdCol', default=NA, help="Name of column with tx id"),
  make_option(c("-a", "--abundanceCol"), type='character', dest='abundanceCol', default=NULL, help="Name of column with abundances (e.g. TPM or FPKM)"),
  make_option(c("-u", "--countsCol"), type='character', dest='countsCol', default=NULL, help="Name of column with estimated counts"),
  make_option(c("-e", "--lengthCol"), type='character', dest='lengthCol', default=NULL, help="Name of column with feature length information"),
  make_option(c("-p", "--importer"), type='character', dest='importer', default=NULL, help="A function used to read in the files"),
  make_option(c("-j", "--existenceOptional"), type='logical', dest='existenceOptional', default=FALSE, help="logical, should tximport not check if files exist before attempting import (default FALSE, meaning files must exist according to file.exists)"),
  make_option(c("-r", "--readLength"), type='numeric', dest='readLength', default=75, help="numeric, the read length used to calculate counts from StringTie's output of coverage. Default value (from StringTie) is 75. The formula used to calculate counts is: cov * transcript length / read length")
)

mandatory <- c("files", "type", "outputCountsFile", "outputAbundancesFile")

opt <- wsc_parse_args(option_list, mandatory = mandatory)

# Read file locations from a file

quantfiles <- readLines(opt$files)

# Name and order the quantification files

names(quantfiles) <- basename(dirname(quantfiles))

# Read transcript to gene mapping

tx2gene <- NA
if ( ! is.null(opt$tx2gene)){
    tx2gene <- read.csv(opt$tx2gene)
}

importer <- opt$importer
if ( ! is.null(importer)){
    importer <- get(importer)
}

# Run tximport

txi <- tximport(quantfiles, type = opt$type, txIn=opt$txIn, txOut=opt$txOut, countsFromAbundance=opt$countsFromAbundance, tx2gene=tx2gene, varReduce=opt$varReduce, dropInfReps=opt$dropInfReps, ignoreTxVersion=opt$ignoreTxVersion, ignoreAfterBar=opt$ignoreAfterBar, geneIdCol=opt$geneIdCol, txIdCol=opt$txIdCol, abundanceCol=opt$abundanceCol, countsCol=opt$countsCol, lengthCol=opt$lengthCol, importer = importer, existenceOptional=opt$existenceOptional, readLength=opt$readLength)

## Generate outputs

# If sparse required, load the requisite libraries

if ( opt$outputFormat != 'tsv' ){
  suppressPackageStartupMessages(require(Matrix))
  suppressPackageStartupMessages(require(DropletUtils))

  write10xCounts(dirname(opt$outputCountsFile), as(txi$counts, "dgCMatrix"), gene.id=rownames(txi$counts), barcodes=colnames(txi$counts), overwrite=TRUE)    
  write10xCounts(dirname(opt$outputAbundancesFile), as(txi$abundance, "dgCMatrix"), gene.id=rownames(txi$counts), barcodes=colnames(txi$abundance), overwrite=TRUE)    
}else{
  write.table(data.frame(Feature=rownames(txi$counts), txi$counts), file=opt$outputCountsFile, sep='\t', quote=FALSE, row.names = FALSE)
  write.table(data.frame(Feature=rownames(txi$abundance), txi$abundance), file=opt$outputAbundancesFile, sep='\t', quote=FALSE, row.names = FALSE)
}

# Read and process the stats files from Kallisto. Equivalent for other methods
# not yet implemented

if ( opt$type == 'kallisto' && ! is.na(opt$outputStatsFile)){
    meta_info_files <- unlist(lapply(quantfiles, function(x) file.path(dirname(x), 'run_info.json')))
    names(meta_info_files) <- names(quantfiles)
    stats <- do.call(rbind, lapply(meta_info_files, function(mif){
        json <- fromJSON(file = mif)
        unlist(lapply(json, function(j) paste(j, collapse = ' ')))
    }))
    write.table(stats, file = opt$outputStatsFile, quote = FALSE, sep="\t")
}
