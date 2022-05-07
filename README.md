# STARR-seq
STARR-seq data analyze

## QC
QC pipeline for MPRA and STARR-seq experiments
> https://github.com/ENCODE-AWG/mpra_starr_qc

- Create fragment count files from STARR-seq BAM files.
  ```R
 #!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(optparse)
    library(data.table)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(RColorBrewer)
    library(doParallel)
    library(argparse)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-i", "--input", nargs="+", help="BAM files with reads")
#parser$add_argument("-p", "--cores", default=4, help="Number of cores to use in multicore processing")
parser$add_argument("--UMI", required=F, action="store_true", default=FALSE, help="Specify if UMI/barcode present for fragments")
parser$add_argument("-o", "--outfile", required=F, help="Output filename to create bed file.")

args <- parser$parse_args()


registerDoParallel(cores=args$cores);

# Shorthand for "pretty number" formatting
pn = function(value) {
    prettyNum(value, big.mark=",")
}

# Shorthand to print the full argument list
msgout = function(...) {
    write(paste(...), stdout());
}

if (!(args$UMI)){
    bfilters = ScanBamParam(mapqFilter=10, flag=scanBamFlag(isSecondaryAlignment=F));
} else {
    bfilters = ScanBamParam(mapqFilter=10, flag=scanBamFlag(isSecondaryAlignment=F, isDuplicate=F));
}

count_reads = function(reads) {
    uniq = unique(reads);
    # sum over duplicates to get a count for each unique 5'/3' end
    uniq$count = countOverlaps( uniq, reads, type="equal" );
    return( uniq );
}

yield.bam = function(X) {
    y = GRanges( readGAlignmentPairs(X, use.names=F, param=bfilters ));
    return(y);
}

map.bam = function(X) {
    return(X);
}

reduce.bam = function(x, y) {
    x = append(x, y);
    # print the number of readpairs processed
    msgout(pn(length(x)), 'mapped human reads');
    return(x);
}

merge_counts = function( x, y, name ) {
    if(length(x)) {
        hits = findOverlaps(x, y, type="equal");
        mcols(x)[hits@from, name] = mcols(x)[hits@from, name] + mcols(y)[hits@to, name];
        y = y[-hits@to,];
        for( cn in colnames(mcols(x))) {
            if( !cn %in% colnames(mcols(y)) ) {
                mcols(y)[,cn] = 0;
            }
        }
        mcols(y) = mcols(y)[,colnames(mcols(x))];
        colnames(mcols(y)) = colnames(mcols(x));
    }
    return(append(x, y));
}


ctFile = args$input;
msgout( "Processing ", ctFile );
infile = BamFile(ctFile, yieldSize=1 * 10^6, asMates=T );
aligned = reduceByYield( infile, yield.bam, map.bam, reduce.bam, parallel=F );

#msgout(pn(length(aligned)), 'mapped reads');

# compute coverage from identical reads => 'count' column
seqlib = count_reads(aligned);
print(seqlib)

# convert the GRanges object to dataframe and add columns corresponding to name, score and barcode information as per the common file format needed

seqlib = as.data.frame(seqlib)
cols_to_add = data.frame(name=paste0(seqlib$seqnames,"_", seqlib$start,"_",seqlib$end), 
                         score = pmin(seqlib$count,1000), 
                         barcode = '.')

seqlib_df = cbind(seqlib, cols_to_add)
col_order = c('seqnames','start','end','name','score','strand','count','barcode')
seqlib_df = seqlib_df[, col_order]

# write the dataframe into output format required

write.table(seqlib_df, file=paste0(args$outfile), quote=FALSE, row.names=FALSE, sep='\t');
  ```
- Aggregate fragments in fixed size bins.
```R
suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(optparse)
    library(data.table)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(RColorBrewer)
    library(doParallel)
    library(argparse)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-i", "--input", nargs="+", help="Input count tables at fragment-level")
#parser$add_argument("-p", "--cores", default=4, help="Number of cores to use in multicore processing")
parser$add_argument("--csaw", required=F, action="store_true", default=FALSE, help="Specify if CSAW package will be used for enhancer call")
parser$add_argument("--starrpeaker", required=F, action='store_true', default=FALSE, help="Specify if STARRPeaker package will be used for enhancer call")
parser$add_argument("--bin", required=F, action='store_true', default=FALSE, help="Specify if fragments need to be binned")
parser$add_argument('-s', "--size", required=F, default=50, type="integer", help="Specify bin size")
#parser$add_argument('--ori', required=F, action="store_true", default=FALSE, help="Specify if orientation needs to be separated")
parser$add_argument("-o", "--outfile", required=T, help="Output filename")
parser$add_argument("--minFragmentLength", required=F, default = -1, type="integer", help="Specify the minimum size of fragments")
parser$add_argument("--maxFragmentLength", required=F, default = -1, type="integer", help="Specify the maximum size of fragments")


# this function will be used to merge counts from replicates
merge_counts = function( x, y, name ) {
    if(length(x)) {
        hits = findOverlaps(x, y, type="equal");
        mcols(x)[hits@from, name] = mcols(x)[hits@from, name] + mcols(y)[hits@to, name];
        y = y[-hits@to,];
        for( cn in colnames(mcols(x))) {
            if( !cn %in% colnames(mcols(y)) ) {
                mcols(y)[,cn] = 0;
            }
        }
        mcols(y) = mcols(y)[,colnames(mcols(x))];
        colnames(mcols(y)) = colnames(mcols(x));
    }
    return(append(x, y));
}

# this function will merge (add) metadata columns for
# identical reads. used to bin data (if desired)
collapse_reads = function( x ) {
  u = x[!duplicated(x)];
  x = x[ duplicated(x)];
  gc();
  for( cn in 1:ncol(mcols(u)) ) {
    z = x[ mcols(x)[,cn]>0 ];
    hits = findOverlaps( u, z, type="equal" );
    if( any(duplicated(hits@from)) ) {
        hits = aggregate( mcols(z)[hits@to,cn] ~ hits@from, FUN='sum' );
        mcols(u)[hits[,1],cn] = mcols(u)[hits[,1],cn] + hits[,2];
    } else {
        mcols(u)[hits@from,cn] = mcols(u)[hits@from,cn] + mcols(z)[hits@to,cn];
    }
  }
  return(u);
}

bin_fragments = function(count, binSize=50){
    
    start(count) = binSize * as.integer(start(count)/ binSize)+1
    end(count) = binSize * as.integer(end(count)/ binSize)+1
    count = collapse_reads(count)
    return(count)
    
}

args <- parser$parse_args()

# Merging count from replicates 
frag_count = GRanges()

for( i in 1:length(args$input) ){
    x = read.csv(args$input[i], sep='\t')
    x = x[, c('seqnames', 'start', 'end', 'count')]
    
    # set sample name
    sample = strsplit(args$input[i], "/")[[1]]
    sample = strsplit(sample[length(sample)], ".csv")[[1]]
    
    # Convert dataframe to GRanges
    x = makeGRangesFromDataFrame(x, 
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,seqnames.field=c("seqnames", "seqname","chromosome", "chrom","Chr", "chromosome_name","seqid"),
                               start.field="Start",
                               end.field=c("End", "stop"),
                               strand.field="Strand",
                               starts.in.df.are.0based=FALSE) # should they?
    colnames(mcols(x)) = c(sample)
    
    # Filter fragments based on size
    if(args$minFragmentLength != -1 & args$maxFragmentLength != -1){
        x = x[width(x) >= args$minFragmentLength & width(x) <= args$maxFragmentLength]
    }
    
    mcols(frag_count)[,sample] = 0;
    frag_count = merge_counts(frag_count, x, sample);

}

# Bin fragments

if(args$bin){
    binned_frag_count = bin_fragments(frag_count, binSize=args$size)
#     # Save
#     write.table(binned_frag_count, file=paste0(args$outfile), quote=FALSE, sep='\t');
    
    # Convert binned fragments count GRanges to dataframe for each replicate 
    for( i in 1:length(colnames(mcols(binned_frag_count)))){
        sample = colnames(mcols(binned_frag_count))[i]
        seqlib = mcols(binned_frag_count)[, sample]
        tmp = binned_frag_count
        mcols(tmp) = seqlib
        seqlib = as.data.frame(tmp) 
        seqlib = seqlib[, c(1,2,3,5,6)]
        colnames(seqlib) = c('chr','start','end', 'strand', 'count')

        cols_to_add = data.frame(name=paste0(seqlib$chr,"_", seqlib$start,"_",seqlib$end), 
                         score = '.', 
                         barcode = '.')

        seqlib_df = cbind(seqlib, cols_to_add)
        col_order = c('chr','start','end','name','score','strand','count','barcode')
        seqlib_df = seqlib_df[, col_order]
        
        # Save
        write.table(seqlib_df, file=paste0(args$outfile[i]), quote=FALSE, sep='\t');
    }
}
```
- Main script to generate QC plots and summary statistics.
```R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(argparse)
})

main_wrap = function(raw_DNA_folder, raw_RNA_folder, threshold.reads=-Inf, 
                     method="MPRA", output.name="qc_output") {
  
  # TO-DO: Add STARR option
  if(method!="MPRA") stop("No current implementation for non-MPRA experiments!")
  
  ### glossary:
  # read - numbers in column 7, usually number of raw sequencing reads
  # unit - column 8 (if provided, currently ignored), UMIs/regions/barcodes
  # element - column 4 (aka sequence, oligo, insert, tile...)
  
  # find files
  raw_DNA_files = list.files(raw_DNA_folder, pattern = '.bed(.gz)?',full.names = T)
  raw_RNA_files = list.files(raw_RNA_folder, pattern = '.bed(.gz)?',full.names = T)
  
  ### check number of samples
  DNA_len = length(raw_DNA_files)
  RNA_len = length(raw_RNA_files)
  if(RNA_len==0 | DNA_len ==0) stop("Couldn't find DNA or RNA files (BED or BED.GZ)")
  if(RNA_len==DNA_len) {
    message("Found ", DNA_len, " DNA and RNA replicates.") } else {
      if(RNA_len>DNA_len & DNA_len==1) {
        message("Found only one DNA reference for ", RNA_len, " RNA replicates.") } else {
          stop("Found ", list.files(raw_DNA_folder, pattern = '.bed(.gz)?'),
               " DNA samples and ", list.files(raw_RNA_folder, pattern = '.bed(.gz)?'), 
               " RNA samples. That's unacceptable, aborting.")
        }
    }
  
  ### make output folder
  output.dir = file.path("report",output.name)
  dir.create(output.dir,showWarnings = F,recursive = T)
  # round(file.size(raw_DNA_files[replicate])/1024^2,1) # file size in MB, if needed
  
  ### extract stats
  stats = list(); ratios=list()
  recalculate_threshold = is.infinite(threshold.reads)
  for(replicate in seq_len(RNA_len)) {
    
    message("Replicate ",replicate)
    
    # for special case of DNA_len==1 and RNA_len>1, do not read DNA multiple times
    if(!(DNA_len==1 & exists("already_read_dna"))) {
      
      if(method == "STARR") {
        # unwritten function to bin STARR regions
        # ideally, should only change the 'name' (4) column
        # consider binning without reading whole table into memory
        STARR_handler(raw_DNA_files[replicate])
        DNA_file = name_of_STARR_handler_file
      } else {
        # for MPRA just read in the raw_mpra format
        DNA_file = raw_DNA_files[replicate]
      }
      
      # consider data.table::fread for large files
      dna = read_tsv(DNA_file,
                     col_names = c('chr','start','end','name','score','strand',
                                   'DNA','barcode'), 
                     # skip header, if it's there
                     skip = as.integer(grepl("start|end|strand",
                                             readLines(raw_DNA_files[replicate],n = 1))),
                     col_types = "ciicicic") %>% 
        group_by(name) %>% 
        summarise(DNA.reads = sum(DNA),
                  DNA.units = n())
      
      ### thresholding: currently only total DNA reads, not units (ie UMI, regions etc)
      # if no threshold provided, use Q1-1.5*IQR
      if(recalculate_threshold) {
        threshold.reads = as.integer(10^(quantile(log10(dna$DNA.reads), probs = c(0.25)) - 
                                     1.5*IQR(log10(dna$DNA.reads))))
      }
      dna_filter = dna %>% filter(DNA.reads>=threshold.reads)
      already_read_dna = T
      
      ### histogram of DNA reads per element
      hist = dna %>%
        select(DNA.reads,name) %>% # units currently ignored
        gather('key','value',-name) %>% # rewrite as pivot_longer
        ggplot(aes(value))+
        geom_histogram(bins=50) +
        geom_vline(data = data.frame(xintercept = c(threshold.reads),
                                     key = c("DNA.reads")), 
                    aes(xintercept=xintercept), color="red") +
        facet_grid(~key) +
        scale_x_log10() +
        labs(title = paste0("DNA replicate ", replicate))+
        xlab('DNA reads or units per element')
      ggsave(file.path(output.dir, paste0("DNA_rep",replicate,"_hist_plot.png")),
             hist, w=4,h=2)
      
    }
    
    rna = read_tsv(raw_RNA_files[replicate],
                   col_names = c('chr','start','end','name','score','strand',
                                 'RNA','barcode'), 
                   # skip header, if it's there
                   skip = as.integer(grepl("start|end|strand",
                                           readLines(raw_RNA_files[replicate],n = 1))),
                   col_types = "ciicicic") %>% 
      group_by(name) %>% 
      summarise(RNA.reads = sum(RNA),
                RNA.units =n())
    
    ### calculate ratios (assuming "DNA is right", ie left_join)
    ratios[[replicate]] = dna_filter %>%
      left_join(rna,by="name") %>%
      mutate(DNA.norm = log2( (DNA.reads+1)/sum(DNA.reads)*1e6 ),
             RNA.norm = log2( (RNA.reads+1)/sum(RNA.reads)*1e6 ),
             log.ratio = RNA.norm - DNA.norm,
             replicate)
    
    ### histogram of RNA reads per element
    hist2 = ratios[[replicate]] %>%
      select(RNA.reads, name) %>% 
      replace_na(list(RNA.reads=1)) %>% 
      gather('key','value',-name) %>% 
      ggplot(aes(value)) +
      geom_histogram(bins=50) +
      facet_grid(~key) +
      scale_x_log10() +
      labs(title = paste0("RNA replicate ", replicate))+
      xlab('RNA reads or units per element')
    ggsave(file.path(output.dir, paste0("RNA_rep",replicate,"_hist_plot.png")),
           hist2, w=4,h=2)
    
    stats[[replicate]] = tibble(output.name, 
                                replicate, 
                                dna_read_threshold = threshold.reads,
                                saturation_rna = round(1-sum(rna$RNA.units)/sum(rna$RNA.reads),2), 
                                saturation_dna = round(1-sum(dna$DNA.units)/sum(dna$DNA.reads),2), 
                                dna_elements_total = nrow(dna), # number of elements
                                dna_elements_filtered = nrow(dna_filter), # number of elements with reads above threshold
                                dna_reads_total = sum(dna$DNA.reads), # number of reads
                                dna_reads_filtered = sum(dna_filter$DNA.reads), # numbers of reads in filtered elements
                                dna_units_filtered = sum(dna_filter$DNA.units), # number of units, e.g. UMIs
                                med_units_per_elem_filtered = median(dna_filter$DNA.units))
  }
  
  # write result to files
  write_tsv(bind_rows(stats),file.path(output.dir,"complete_stats.tsv"))
  complete_ratios = bind_rows(ratios) %>% 
    mutate(label = paste0("rep",replicate))
  write_tsv(complete_ratios,gzfile(file.path(output.dir,"complete_ratios.tsv.gz")))
  

  ### make correlation plots (RNA,DNA,ratios)
  # suppressing warnings because of pairwise incomplete observations
  pl.RNA = complete_ratios %>% 
    select(label,RNA.norm,name) %>% 
    spread(label,RNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "RNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_RNA.png"),pl1,w=7,h=7))
  
  pl.DNA = complete_ratios %>% 
    select(label,DNA.norm,name) %>% 
    spread(label,DNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "DNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_DNA.png"),pl2,w=7,h=7))
  
  pl.ratio = complete_ratios %>% 
    select(label,log.ratio,name) %>% 
    spread(label,log.ratio) %>% 
    select(-name) %>% 
    ggpairs(title = "ratio correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_ratio.png"),pl3,w=7,h=7))
  
}  


# Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-d", "--dnas", required=T, help="Folder with Raw DNA counts")
parser$add_argument("-r", "--rnas", required=T, help="Folder with Raw RNA counts")
parser$add_argument("--method", required=F, choices=c("MPRA", "STARR"),  
                    default="MPRA", help="Specify the experiment type.")
parser$add_argument("--threshold-reads", required=F,   
                    default="-Inf", help="Threshold used to filter out low DNA count elements.")
parser$add_argument("-o", "--outfile", required=F, 
                    help="Output rootname for newly generated QC files.")

args <- parser$parse_args()
# load libraries after successfully parsing the command-line arguments
suppressPackageStartupMessages({
    library(tidyverse)
    library(GGally)
})
theme_set(theme_bw())
main_wrap(args$dnas, args$rnas, 
          threshold.reads=as.numeric(args$threshold_reads), 
          method=args$method, 
          args$outfile)
```
