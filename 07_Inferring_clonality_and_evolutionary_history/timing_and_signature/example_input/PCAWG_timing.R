# Mutational timing of gains in single samples

# Libraries and additional code
library(VariantAnnotation)
library(cancerTiming)
source("cancerTiming.WGD.R")

# each VCF file is input as argument to script  
args = commandArgs(TRUE)
vcf_filename = args[1]

# clustering results
clustering_output_path = args[2]

# copy number
cn_path = args[3]

# samplename from VCF
id = gsub(".vcf","",basename(vcf_filename))
print(id)

# function to assign mutations as clonal or subclonal
assignMuts <- function(subclones_file){
  clonal.rows = NULL
  # define clonal cluster and subclones
  if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) == 1){
    
    # if there is only one clone between 0.9 and 1.1
    # take this as the clone, and add the superclonal mutations
    # subclones are everything below 
    clonal.row = subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]
    clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
    subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
    no.subclones = nrow(subclonal.rows)
    
    # if there is nothing in the range of 0.9 to 1.1
    # remove superclone
    # and count only subclonal clusters and mutations 
  } else if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) == 0){
    
    subclones = subset(subclones, fraction_cancer_cells < 1)
    clonal.rows = as.data.frame(matrix(ncol=3, nrow=0))
    subclonal.rows = subset(subclones, fraction_cancer_cells < 1)
    no.subclones = nrow(subclonal.rows)
    
    # if there are multiple clusters in the range of 0.9 and 1.1 
    # look between 0.95 and 1.05
  } else if (nrow(subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]) > 1){
    
    clonal.row.strict = subclones[subclones$fraction_cancer_cells > 0.95 & subclones$fraction_cancer_cells < 1.05,]
    
    # if there is only 1 cluster in the strict range 
    if (nrow(clonal.row.strict) == 1){
      
      clonal.row = clonal.row.strict
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
      
      # if there are still two clusters in the strict range
      # then take the larger ccf 
    } else if (nrow(clonal.row.strict) > 1){
      
      clonal.row = clonal.row.strict[which.max(clonal.row.strict$fraction_cancer_cells),]
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
      
      # if there is nothing between 0.95 and 1.05
      # but two within 0.9 and 1.1
    } else if (nrow(clonal.row.strict) == 0){
      
      clonal.row.relaxed = subclones[subclones$fraction_cancer_cells > 0.9 & subclones$fraction_cancer_cells < 1.1, ]
      clonal.row = clonal.row.relaxed[which.max(clonal.row.relaxed$n_snvs),]
      clonal.rows = subset(subclones, fraction_cancer_cells >= clonal.row$fraction_cancer_cells)
      subclonal.rows = subset(subclones, fraction_cancer_cells < clonal.row$fraction_cancer_cells)
      no.subclones = nrow(subclonal.rows)
    }
  }
  return(clonal.rows)
}


# these are all files available for the clustering
clustering_output_files = list.files(clustering_output_path, pattern = "cluster.xls", recursive=FALSE, full.names = TRUE)
#subclonal_structure_files = list.files(clustering_output_path, pattern = "subclonal_structure.txt.gz", recursive=FALSE, full.names = TRUE)

# the file containing tumour purity and ploidy
purity_ploidy = read.table("purity.ploidy.txt", header = TRUE)
sample_purity = subset(purity_ploidy, samplename == id)$purity
norm_cont = 1 - sample_purity

# the copy number files
cn_files = list.files(cn_path, recursive=FALSE, full.names = TRUE)  

# get mutations and read counts from vcf
sample_vcf = vcf_filename

# first read in vcf file and format
vcf = readVcfAsVRanges(sample_vcf, "hg19",  param = ScanVcfParam(fixed=c("ALT","FILTER"),geno=NA))
#vcf_df=read.table(vcf_filename, header=T)
print(sample_vcf)
vcf_df = as.data.frame(vcf)
vcf_df = subset(vcf_df, !is.na(vcf_df$t_alt_count) & !is.na(vcf_df$t_ref_count))
#vcf_df = vcf_df[,c(1,2,3,6,7)]

# for each sample, get the clustering output 
sample_clustering_output = clustering_output_files[grep(id, clustering_output_files)]
print(sample_clustering_output)

# get DP output files 
if (file.exists(sample_clustering_output)){
  cl_output = read.table(sample_clustering_output, header=TRUE, sep="\t")
  print("1. Sample has clustering output")
  
  # get subclonal structure file 
  
    clonal_mutations = subset(cl_output, mut_type == "SNV")
    clonal_mutations = subset(clonal_mutations, Clonal == "clonal")
    
    if (is.data.frame(clonal_mutations) == TRUE & nrow(clonal_mutations) > 0){
      print("3. Sample has clonal mutations")
      
      # convert vcf and clonal dp_output to granges, subset vcf to get only clonal mutations
      vcf_gr = makeGRangesFromDataFrame(vcf_df, keep.extra.columns = T)
      clonal_mutations$start = clonal_mutations$position
      clonal_mutations$end = clonal_mutations$position
      clonal_mutations$position = NULL
      clonal_mutations_gr = makeGRangesFromDataFrame(clonal_mutations)
      
      # get overlap between vcf and clonal mutations
      hits = findOverlaps(vcf_gr, clonal_mutations_gr, type = "start")
      idx = hits@from
      vcf_clonal_gr = unique(vcf_gr[idx])
      print(c("vcf_clonal_gr", length(idx)))
      
      # read in consensus copy number file and get clonal regions
      sample_cn = cn_files[grep(id, cn_files)]
      print(sample_cn)
      
      if (file.exists(sample_cn)) {
        print("4. Sample has copy number")
        cn = read.table(sample_cn, header=TRUE)
        
        # get clonal gain segments and reduce file 
        # take only a-d segments 
        cn_gains = subset(cn, major_cn > 1)[,1:6]
        
        if (is.data.frame(cn_gains) == TRUE & nrow(cn_gains) > 0){
          print(c("5. Sample has clonal gains", nrow(cn_gains)))
          
          # make GRanges of clonal gain position
          cn_gains_gr = makeGRangesFromDataFrame(cn_gains, keep.extra.columns=TRUE)
          
          # get mutations in regions of clonal copy number 
          vcf_gains_gr = mergeByOverlaps(vcf_clonal_gr, cn_gains_gr, type = "within")
          muts_df = as.data.frame(vcf_gains_gr$vcf_clonal_gr)
          cn_df = as.data.frame(vcf_gains_gr$cn_gains_gr, row.names=NULL)
	  names_cn_df = names(cn_df)
	  names(cn_df)[which(names_cn_df=="start")] = "segment.start"
          names(cn_df)[which(names_cn_df=="end")] = "segment.end"
          names(cn_df)[which(names_cn_df=="width")] = "segment.width"
          names(cn_df)[which(names_cn_df=="strand")] = "segment.strand"
          
          
          muts_clonal_gains = cbind(muts_df, cn_df)
	  print(names(muts_clonal_gains))
         
          if (is.data.frame(muts_clonal_gains) == TRUE & nrow(muts_clonal_gains) > 0){
            print(c("6. Sample has clonal mutations in clonal gains", nrow(muts_clonal_gains)))
            
            muts_clonal_gains$type = "none"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 2 & muts_clonal_gains$minor_cn == 1] = "SingleGain"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 2 & muts_clonal_gains$minor_cn == 0] = "CNLOH"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 3 & muts_clonal_gains$minor_cn == 1] = "DoubleGain"
            muts_clonal_gains$type[muts_clonal_gains$major_cn == 2 & muts_clonal_gains$minor_cn == 2] = "WGD"
	    print(names(muts_clonal_gains))
            names(muts_clonal_gains)[1] = "chr"
            
            muts_clonal_gains$segId = paste0(muts_clonal_gains$chr,"_",muts_clonal_gains$segment.start, "_", muts_clonal_gains$segment.end, "_", muts_clonal_gains$type)
            muts_clonal_gains = subset(muts_clonal_gains, type != "none")
            
            if (is.data.frame(muts_clonal_gains) == TRUE & nrow(muts_clonal_gains) > 0){
              print(c("7. Sample has clonal mutations within clonal timeable gains", nrow(muts_clonal_gains)))
              muts_clonal_gains$sample = id
              
              # do formatting for eventTiming  
	      muts_clonal_gains$nMutAllele = muts_clonal_gains$t_alt_count
              muts_clonal_gains$nReads = muts_clonal_gains$t_ref_count + muts_clonal_gains$t_alt_count
              muts_clonal_gains$mutationId = paste0(muts_clonal_gains$chr, "_", muts_clonal_gains$start)  
              write.table(muts_clonal_gains, file=paste0("out/", id, "_mcg.txt"), sep="\t", quote = FALSE, row.names=FALSE)
              clonal_events_list = split(muts_clonal_gains, muts_clonal_gains$sample)
              
              arguments = list(bootstrapCI="nonparametric", minMutations=2)
              
              x = eventTimingOverList.WGD(dfList = clonal_events_list, normCont = norm_cont, eventArgs = arguments)
              
              # format results
              y = getPi0Summary(x, CI=TRUE)
              piSum = na.omit(y)
              rownames(piSum) = NULL
              
              print("8. Writing output")
              
              # save as gz files
              write.table(piSum, file=paste0("out/", id, "_timed_segments.txt"), sep="\t", quote = FALSE, row.names=FALSE) 
              
              
            }
          }
        }
      }
    }
} 


