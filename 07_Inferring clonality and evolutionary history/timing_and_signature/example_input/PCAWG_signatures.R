## PCAWG-11 Mutational Signatures timing

## Libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(nnls)
library(reshape2)
library(plyr)
library(scales)

## Input files

# Sample ID
args = commandArgs(TRUE)
id = args[1]
print(id)

# Mutation clustering results
clustering_file_path = args[2]

# VCF files
vcf_file_path = args[3]
indel_file_path = args[4]


## Functions

# Determines if individual SNVs are part of an MNV
isMNV <- function(vcf) {
  d <- diff(start(vcf)) == 1 # are they next to each other?
  w <- c(FALSE, d) | c(d, FALSE) # need to add comparisons at the end, set to FALSE
  return(w)
}

# Formatting of SNVs -> pyrimidine context
formatSNVs <- function(vcf){
  vcf=vcf
  vcf = vcf[!is.na(vcf$t_alt_count) & !is.na(vcf$t_ref_count)] # remove variants with NA for alt or ref
  k=3
  ranges = resize(vcf, k, fix = "center") # resize to get 1 bp either side 
  seqlevelsStyle(ranges) = "UCSC"
  ref=FaFile("/hg19/hg19AddVirus.fa")
  context = getSeq(ref, ranges) # get 3 bp contexts
  ct = as.data.frame(context)
  vcf$context = ct$x
  
  rev.context = reverseComplement(DNAStringSet(ct$x)) # get reverse complement of context as well
  vcf$rev.context = strsplit(toString(rev.context), ", ")[[1]]
  
  # formatting of bp change and contexts
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "A") | (ref(vcf) == "G" & alt(vcf) == "T")] = "C>A"
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "G") | (ref(vcf) == "G" & alt(vcf) == "C")] = "C>G"
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "T") | (ref(vcf) == "G" & alt(vcf) == "A")] = "C>T"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "A") | (ref(vcf) == "A" & alt(vcf) == "T")] = "T>A"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "C") | (ref(vcf) == "A" & alt(vcf) == "G")] = "T>C"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "G") | (ref(vcf) == "A" & alt(vcf) == "C")] = "T>G"
  
  pyrimidine = c("C", "T")
  
  vcf$base_before = ifelse(ref(vcf) %in% pyrimidine, substring(vcf$context, 1, 1), substring(vcf$rev.context, 1, 1))
  vcf$base_after = ifelse(ref(vcf) %in% pyrimidine, substring(vcf$context, 3), substring(vcf$rev.context, 3))
  
  vcf$tnc = paste0(vcf$base_before,"[",vcf$change,"]",vcf$base_after)
  return(vcf)
}

# Format MNVs so that they match the signatures annotation
formatMNVs <- function(vcf){

  u <- reduce(granges(vcf[which(isMNV(vcf)),]))
  u <- u[width(u)==2]
  
  if(length(u)==0) return(u)
  
  seqlevelsStyle(u) = "UCSC"
  
  # Check timing is the same
  t = read.delim(clustering_file_path, 
                 header=TRUE, stringsAsFactors = FALSE)
  t$start = t$position
  t$end = t$position
  t = subset(t, mut_type != "SV" & mut_type != "indel")
  t = makeGRangesFromDataFrame(t, keep.extra.columns = TRUE)
  seqlevelsStyle(t) = "UCSC"
  
  # Merge VCF with timing annotation
  timed_mnvs = mergeByOverlaps(t, u)
  
  # Remove single variants
  timed_mnvs$lengths = rep(rle(start(timed_mnvs$u))$lengths, rle(start(timed_mnvs$u))$lengths)
  timed_mnvs = subset(timed_mnvs, lengths>1)
  
  if(nrow(timed_mnvs)==0) return(timed_mnvs$u)
  
  times = lapply(split(timed_mnvs$timing, sort(rep(1:(nrow(timed_mnvs)/2), 2))), unique)
  times = !lapply(times, length) > 1
  
  u = reduce(timed_mnvs$u)[times]
  
  if(length(u)==0) return(u)

  seqlevelsStyle(vcf) = "UCSC"
  r <- as.character(ref(vcf))[vcf %over% u] # get reference base
  h <- subjectHits(findOverlaps(granges(vcf),u))
  a <- as.character(alt(vcf))[vcf %over% u] # get alternate base
  rr <- sapply(split(r, h), paste, collapse="")
  aa <- sapply(split(a, h), paste, collapse="")
  
  u$ref_mnv = rr
  u$alt_mnv = aa
  
  convert_all = c("AA","AG","CA","GA","GG","GT")
  
  # change purine bases, this is so they match the bases given in DBS signatures file
  c = which(u$ref_mnv %in% convert_all)
  u$ref_mnv[c] = reverseComplement(DNAStringSet(u$ref_mnv[c])) # if there are no purine bases it just won't add anything
  u$alt_mnv[c] = reverseComplement(DNAStringSet(u$alt_mnv[c]))
  
  # Format
  u$DBS = paste0(u$ref_mnv, ">", u$alt_mnv)
  
  # There are specific cases that need to be converted
  u$DBS[u$DBS=="AT>GG"] = "AT>CC"
  u$DBS[u$DBS=="AT>TC"] = "AT>GA"
  u$DBS[u$DBS=="AT>TG"] = "AT>CA"
  u$DBS[u$DBS=="CG>AA"] = "CG>TT"
  u$DBS[u$DBS=="CG>AC"] = "CG>GT"
  u$DBS[u$DBS=="CG>GA"] = "CG>TC"
  u$DBS[u$DBS=="GC>CT"] = "GC>AG"
  u$DBS[u$DBS=="GC>TG"] = "GC>CA"
  u$DBS[u$DBS=="GC>TT"] = "GC>AA"
  u$DBS[u$DBS=="TA>AC"] = "TA>GT"
  u$DBS[u$DBS=="TA>AG"] = "TA>CT"
  u$DBS[u$DBS=="TA>CC"] = "TA>GG"
  
  return(u)
}

seqDiff <- function(x,y,diff){
  x1 = unlist(strsplit(x, ""))
  y1 = unlist(strsplit(y, "")) 
  if (diff == "deleted"){
    return(x1[!1:length(x1) %in% pmatch(y1, x1)])
  } 
  if (diff == "inserted"){
    return(y1[!1:length(y1) %in% pmatch(x1, y1)])
  }
}


# Get indels that correspond to characteristic events from indel signatures ID1, ID2, ID13 and ID8
formatIndels <- function(vcf){
  
  indel_vcf = vcf
  
  if(length(indel_vcf)<=1) return(indel_vcf)
  
  seqlevelsStyle(indel_vcf) = "UCSC"
  
  # Start by computing what exactly is changed in the indel
  # NB - Some are both insertions and deletions e.g. AAA -> TT
  
  # Are all the reference bases in the alternate allele, in the right order?
  indel_vcf$deleted = !mapply(grepl, ref(indel_vcf), alt(indel_vcf))
  indel_vcf$inserted = !mapply(grepl, alt(indel_vcf), ref(indel_vcf))
  indel_vcf$del_base = mapply(seqDiff, ref(indel_vcf), alt(indel_vcf), diff="deleted", SIMPLIFY = FALSE)
  indel_vcf$ins_base = mapply(seqDiff, ref(indel_vcf), alt(indel_vcf), diff="inserted", SIMPLIFY = FALSE)
  
  # We are looking for 4 different types of indels, corresponding to ID1, ID2, ID8 and ID13
  
  # 1. ID1: 1bp insertion of T at regions of 5+ T's, or A at region of 5+ A's
  ins_1bp_T = indel_vcf[indel_vcf$ins_base=="T" | indel_vcf$ins_base=="A"]
  
  if (length(ins_1bp_T)>0){
    
    # Get 11bp surrounding, rle of T's/ A's, are they longer or equal to 5? 
    k=11
    ins_1bp_T = resize(ins_1bp_T, k, fix = "center") # get 5 bp either side of reference
    ref=FaFile("/hg19/hg19AddVirus.fa")
    context_ins_T = as.data.frame(getSeq(ref, ins_1bp_T))
    context_ins_T$x = substring(context_ins_T$x, 2) # take the first base off, so that 5bp sequences must include centre base or be on the right hand side
    
    ins_1bp_T_5hp = c()
    for (n in 1:length(ins_1bp_T)){
      r = rle(strsplit(context_ins_T$x[n], "")[[1]]) # go through and count how many consecutive bases
      if (any(r$lengths >= 5)){ # are there any stretches >= 5?
        b = r$values[which(r$lengths>=5)] # get the 5+ repeated base
        ins_b = ins_1bp_T$ins_base[n] # and the inserted base
        if (ins_b %in% b){ # is inserted bp one of the repeated ones? In this case, this is fine to be either
          ins_1bp_T_5hp = c(ins_1bp_T_5hp, n) # if they are (should be A and A or T and T, then add this number to vector)
        }
      } 
    }
    
    ID1 = ins_1bp_T[ins_1bp_T_5hp]
    
    if (length(ID1)>0){ 
      ID1$change = NULL
      k=1
      ID1 = resize(ID1, k, fix="center") # resize back to centre single position
      ID1$classification = "ID1" 
    } # otherwise it is just the empty granges to be added
    
  } else {
    
    ID1 = ins_1bp_T # this will be empty, can just add to the others at the end
    
  }
  
  # 2. ID2: 1bp deletion of T at regions of 5+ T's, or A at regions of 5+ A's
  del_1bp_T = indel_vcf[indel_vcf$del_base=="T" | indel_vcf$del_base=="A"]
  
  if (length(del_1bp_T)>0){
    
    # Get 11bp surrounding, rle of T's/ A's, are they longer or equal to 6? 
    k=12
    del_1bp_T = resize(del_1bp_T, k, fix = "center") 
    ref=FaFile("/hg19/hg19AddVirus.fa")
    context_del_T = as.data.frame(getSeq(ref, del_1bp_T))
    context_del_T$x = substring(context_del_T$x, 2)
    
    del_1bp_T_6hp = c()
    for (n in 1:length(del_1bp_T)){
      r = rle(strsplit(context_del_T$x[n], "")[[1]]) # go through and count how many consecutive bases
      if (any(r$lengths >= 6)){ # are there any stretches >= 6?
        b = r$values[which(r$lengths>=6)] # get the 6+ repeated base
        if (length(b)>1) { b = b[2] }
        del_b = del_1bp_T$del_base[n] # and the deleted base
        if (del_b %in% b){ # are they the same?
          del_1bp_T_6hp = c(del_1bp_T_6hp, n) # if they are (should be A and A or T and T, then add this number to vector)
        }
      } 
    }
    
    ID2 = del_1bp_T[del_1bp_T_6hp]
    
    if (length(ID2)>0) {
      ID2$change = NULL
      k = 1
      ID2 = resize(ID2, k, fix="center")
      ID2$classification = "ID2"
    }
    
  } else {
    ID2 = del_1bp_T
  }
  
  
  total_IDs = c(ID1, ID2)
  return(total_IDs)
  
} 




# Get timed multinomials of SNVs, MNVs and indels, removing MNVs from SNV calls
getSNVs_MNVs_IDs <- function(id){
  
  # Check in both snv_mnv and graylist subfolders, read in VCF of SNVs
 vcf = readVcfAsVRanges(vcf_file_path, "hg19", param = ScanVcfParam(fixed=c("ALT","FILTER"),geno=NA))
  
  # Get MNVs from the SNV vcf
  mnv_vcf = vcf[which(isMNV(vcf))]
  
  # Remove MNVs from SNVs
  snv_vcf = subsetByOverlaps(vcf, mnv_vcf, invert=TRUE)
  
  # Now get indels
  indel_vcf = readVcfAsVRanges(indel_file_path, "hg19", param = ScanVcfParam(fixed=c("ALT","FILTER"),geno=NA))
  
  # Format each type of event to agree with classification in signatures
  
  # SNVs
  snv_vcf = formatSNVs(vcf=snv_vcf)
  snv = GRanges(seqnames = seqnames(snv_vcf), ranges = ranges(snv_vcf), class=snv_vcf$tnc) # Simplify VCF to just position and mutation class
  seqlevelsStyle(snv) = "UCSC"
  snv$type = "SNV"
  
  # MNVs
  mnv_vcf = formatMNVs(vcf=vcf) # Give it the unfiltered SNV VCF, as MNVs are identified within function
  seqlevelsStyle(mnv_vcf) = "UCSC"
  if (length(mnv_vcf)>0){ # Otherwise just leave mnv_vcf as empty GRanges, will be added to total and contribute nothing
    mnv_vcf$ref_mnv = NULL
    mnv_vcf$alt_mnv = NULL
    mnv_vcf$class = mnv_vcf$DBS
    mnv_vcf$DBS = NULL
    mnv_vcf$type = "SNV"
  }  

  # Indels
  indel = formatIndels(vcf=indel_vcf)
  if(length(indel)>1){
    seqlevelsStyle(indel) = "UCSC"
    indel = GRanges(seqnames = seqnames(indel), ranges = ranges(indel), class=indel$classification)
    rm(indel_vcf)
    indel$type = "indel"
  } else {
    indel = GRanges()
  }

  # Combine
  total_vcf = c(snv, mnv_vcf, indel)
  rm(vcf)
  
  # Get timing annotation
  t = read.delim(clustering_file_path, 
                 header=TRUE, stringsAsFactors = FALSE)
  t$start = t$position
  t$end = t$position
  t = subset(t, mut_type != "SV")
  t = makeGRangesFromDataFrame(t, keep.extra.columns = TRUE)
  seqlevelsStyle(t) = "UCSC"
  
  # Merge VCF with timing annotation
  timed = mergeByOverlaps(t, total_vcf, type="start")
  timed = subset(timed, mut_type==type)
  
  # Double "clonal" section and add
  if ("clonal" %in% timed$timing){
    clonal = subset(timed, timing!="subclonal")
    clonal$timing = "clonal"
    timed = rbind(timed, clonal)
    # Also add total
    all = subset(timed, timing!="clonal")
    all$timing = "total"
    timed = rbind(timed, all)
  } else{
    timed = timed
    all = timed
    all$timing = "total"
    timed = rbind(timed, all)
  }
  
  # Split by timing
  timed = split(timed, timed$timing)
  
  # Make multinomial comprising all classes, per time frame
  timed_multi = lapply(timed, function(x){
    table(x$class)
  })
  
  
}

# Applies NNLS
nnls_sol <- function(m, s_active, no_active) {
  
  a = matrix(0, nrow = 1, ncol = no_active)
  for (i in 1:ncol(m))
    a[i,] = coef(nnls(s_active, m[,i]))
  a  
}

# Get signature weights from mutation matrix
getSignatureWeights <- function(mutationMatrix, active_comp, sig_comp, sigs){
  
  m = as.matrix(mutationMatrix)
  
  if (is.matrix(m) & is.numeric(m[1]) & is.matrix(active_comp) & is.numeric(active_comp[1]) & nrow(active_comp) == nrow(m)){
    
    no_active = ncol(active_comp)
    a = nnls_sol(m,active_comp,no_active)
    
    colnames(a) = sigs
    
    # add missing signatures set to 0 
    a_df = as.data.frame(a)
    add = subset(colnames(sig_comp), !colnames(sig_comp) %in% colnames(a_df))
    vals = rep(0,length(add))
    names(vals) = add
    vals = t(as.data.frame(vals))
    rownames(vals) = NULL
    a_final = cbind(a, vals)
    a_final = as.data.frame(a_final)
    a_final = a_final[,colnames(sig_comp)]
    a_final$sample = id
    
    return(a_final)
    
  } 
}


extractSigs <- function(events){

  # start with SNV signatures
  snv_id = subset(snv_samples, Sample.Name==id)
  active_sigs = names(snv_samples[,2:ncol(snv_samples)][which(subset(snv_samples, Sample.Name==id)[,2:ncol(snv_samples)] > 0)])
  active_comp = as.matrix(snv_comp[,active_sigs])
  events_snv = lapply(events, function(x){
    x = x[names(x) %in% rownames(snv_comp)]
    if (length(names(x))!=96){
      missing = rownames(snv_comp)[!rownames(snv_comp) %in% names(x)]
      add = rep(0, length(missing))
      names(add) = missing
      x = c(x, add)
      x = x[rownames(snv_comp)]
    }
    x = x[rownames(snv_comp)]
    return(x)
  })
  weightsList = lapply(events_snv, getSignatureWeights, active_comp=active_comp, sig_comp=snv_comp, sigs=active_sigs)
  weights_snv = do.call("rbind", weightsList)
  
  # MNV signatures
  mnv_id = subset(mnv_samples, Sample.Name==id) 
  tailn = length(colnames(mnv_id))
  if(nrow(mnv_id)>0 && sum(mnv_id[,2:tailn])>0){ # if the sample has MNV signatures
    active_mnv = names(mnv_samples[,2:tailn][which(subset(mnv_samples, Sample.Name==id)[,2:tailn] > 0)])
    active_mnv_comp =  as.matrix(mnv_comp[,active_mnv])
    events_mnv = lapply(events, function(x){
      x = x[names(x) %in% rownames(mnv_comp)]
      if (length(names(x))!=78){
        missing = rownames(mnv_comp)[!rownames(mnv_comp) %in% names(x)]
        add = rep(0, length(missing))
        names(add) = missing
        x = c(x, add)
        x = x[rownames(mnv_comp)]
      } 
      x = x[rownames(mnv_comp)]
      return(x)
    })
    weightsList_MNV = lapply(events_mnv, getSignatureWeights,active=active_mnv_comp,sig_comp=mnv_comp,sig=active_mnv)
    weights_mnv = do.call("rbind", weightsList_MNV)
  } else {
    weights_mnv = data.frame(matrix(0, nrow=nrow(weights_snv), ncol=12))
    colnames(weights_mnv) = colnames(mnv_id)[2:ncol(mnv_id)]
  }
  
  # ID signatures
  events_id = lapply(events, function(x){
    x = x[grep("ID", names(x))]
    indel_sigs = c("ID1", "ID2")
    missing = indel_sigs[which(!indel_sigs %in% names(x))]
    add = rep(0, length(missing))
    names(add) = missing
    x = c(x, add)
    x = x[indel_sigs]
  })
  events_id = do.call("rbind", events_id)
  if(ncol(events_id)!=2){
    add = data.frame(matrix(0, nrow=nrow(events_id), ncol=length(missing)))
    colnames(add) = missing
    events_id = cbind(events_id, add)
  }
  sample_id = subset(id_comp, Sample.Name==id)
  z = names(sample_id[colnames(sample_id) %in% colnames(events_id)][which(sample_id[colnames(sample_id) %in% colnames(events_id)] == 0)])
  events_id[,colnames(events_id) %in% z] = 0
  
  
  # Combine signature results
  all_weights = cbind(weights_snv[,1:ncol(weights_snv)-1],weights_mnv[,1:ncol(weights_mnv)-1], events_id)
  all_weights$sample = id
  return(all_weights)
}

# Bootstrapping function
replicateSignatureChanges <- function(multi_list, time){
  
  # Events list should have all events, even if 0
  all_events = c(rownames(snv_comp), rownames(mnv_comp), "ID1", "ID2")
  multi_list = lapply(multi_list, function(x){
    add = rep(0, length(all_events[which(!all_events %in% names(x))]))
    names(add) = all_events[which(!all_events %in% names(x))]
    x = c(x, add)
    x[names(x) %in% c("ID1", "ID2")] = x[names(x) %in% c("ID1", "ID2")] + 0.001 
    x
  })
  
  a_List_resample = lapply(multi_list, function(x){ 
    rs = t(rmultinom(1, size=sum(x), prob=x)) 
    x = as.numeric(rs)
    names(x) = colnames(rs)
    x
  })
  
  # get weights from new mutations 
  weights_resample = extractSigs(events=a_List_resample)
  
  
  # make proportional, using same totals as for actual changes
  weights_resample$time_total = timed_muts[match(rownames(weights_resample), names(timed_muts))]
  weights_resample$n_unassigned = weights_resample$time_total - rowSums(weights_resample[,1:totalsign])
  weights_resample[,1:totalsign][,which(colSums(weights_resample[,1:totalsign]) == 0)] = NA

  weights_resample[,c(1:totalsign, totalsign+3)] = weights_resample[,c(1:totalsign, totalsign+3)]/rowSums(weights_resample[,c(1:totalsign, totalsign+3)], na.rm=TRUE)
  weights_resample$time = rownames(weights_resample)
  
  weights_resample$n_unassigned[weights_resample$n_unassigned < 0] = 0
  weights_resample[rowSums(weights_resample[,1:totalsign], na.rm=TRUE) > 1,c(1:totalsign,totalsign+3)] = weights_resample[rowSums(weights_resample[,1:totalsign], na.rm=TRUE) > 1,c(1:totalsign,totalsign+3)]/rowSums(weights_resample[rowSums(weights_resample[,1:totalsign], na.rm=TRUE) > 1,c(1:totalsign,totalsign+3)], na.rm=TRUE)
  
  # Calculate changes from early-late and clonal-subclonal
  all_weights_long_resample = reshape(suppressMessages(melt(weights_resample[,c(1:(totalsign+1), (totalsign+3):(totalsign+4))])),
                                      timevar="time",
                                      idvar=c("sample","variable"),
                                      direction="wide")
  colnames(all_weights_long_resample)[3:ncol(all_weights_long_resample)] = weights_resample$time
  
  # Add pseudocounts and rescale
  all_weights_long_resample[,3:ncol(all_weights_long_resample)] = all_weights_long_resample[,3:ncol(all_weights_long_resample)] + 0.001
  all_weights_long_resample[,3:ncol(all_weights_long_resample)] = sweep(all_weights_long_resample[,3:ncol(all_weights_long_resample)],2,colSums(all_weights_long_resample[,3:ncol(all_weights_long_resample)], na.rm=TRUE),`/`)
  
  
  if (time=="el"){
    # If there are early and late columns, calculate early and late changes
    all_weights_long_resample$early_late = log2((all_weights_long_resample$`clonal [late]`/(1 - all_weights_long_resample$`clonal [late]`)) / (all_weights_long_resample$`clonal [early]`/(1 - all_weights_long_resample$`clonal [early]`)))
    changes = all_weights_long_resample$early_late
    names(changes) = all_weights_long_resample$variable
    changes = as.data.frame.matrix(t(as.data.frame(changes)))
  } else if (time=="cs"){
    # Same for clonal and subclonal
    all_weights_long_resample$clonal_subclonal = log2((all_weights_long_resample$subclonal/(1 - all_weights_long_resample$subclonal)) / (all_weights_long_resample$clonal/(1 - all_weights_long_resample$clonal)))
    changes = all_weights_long_resample$clonal_subclonal
    names(changes) = all_weights_long_resample$variable
    changes = as.data.frame.matrix(t(as.data.frame(changes)))
  }
  
  return(changes)
}




# Read in signatures
snv_comp = read.table("SBS.sig.txt", header = TRUE)
mnv_comp = read.table("DBS.sig.txt", header = TRUE)

# Change formatting
snv_comp$type = snv_comp$Mutation.Type
rownames(snv_comp) = snv_comp$type
snv_comp = snv_comp[,2:length(colnames(snv_comp))]

rownames(mnv_comp) = mnv_comp$Mutation.Type
mnv_comp$Mutation.Type = NULL
MNV_LEVELS = rownames(mnv_comp)

# And their weights in the samples
snv_samples = read.table("SBS.act.txt", header = TRUE)
mnv_samples = read.table("DBS.act.txt", header = TRUE)
id_comp = read.table("ID.act.txt", header = TRUE)

# Get the multinomials
events = getSNVs_MNVs_IDs(id = id)

# Extract the signatures
n_weights = extractSigs(events=events)

# Combine the signatures here

# Add a column of the number of mutations per time frame
timed_muts = read.delim(clustering_file_path, 
                        header=TRUE, stringsAsFactors = FALSE)
#timed_muts$timing = timed_muts$Clonal   ########trans el or cs
timed_muts = table(timed_muts$timing)
timed_muts = c(timed_muts, clonal=sum(timed_muts[grep("clonal ", names(timed_muts))]), total=sum(timed_muts))
totalsign = length(colnames(n_weights))-1
n_weights$time_total = timed_muts[match(rownames(n_weights), names(timed_muts))]
n_weights$n_unassigned = n_weights$time_total - rowSums(n_weights[,1:totalsign])

# Make proportional
all_weights = n_weights

# Where all values are 0, make NA
all_weights[,1:totalsign][,which(colSums(all_weights[,1:totalsign]) == 0)] = NA
all_weights$time = rownames(all_weights)

# Make proportional
all_weights[,c(1:totalsign, totalsign+3)] = all_weights[,c(1:totalsign, totalsign+3)]/all_weights$time_total

# If more mutations are assigned than are actually in the sample, rescale so sum to 1
all_weights$n_unassigned[all_weights$n_unassigned < 0] = 0
all_weights[rowSums(all_weights[,1:totalsign], na.rm=TRUE) > 1, c(1:totalsign,totalsign+3)] = all_weights[rowSums(all_weights[,1:totalsign], na.rm=TRUE) > 1, c(1:totalsign,totalsign+3)]/rowSums(all_weights[rowSums(all_weights[,1:totalsign], na.rm=TRUE) > 1, c(1:totalsign,totalsign+3)], na.rm=TRUE)

write.table(all_weights, paste0(id, "_all_weights.txt"), sep="\t", quote=FALSE)

# Calculate changes from early-late and clonal-subclonal
all_weights_long = reshape(melt(all_weights[,c(1:(totalsign+1), (totalsign+3):(totalsign+4))]),
                           timevar="time",
                           idvar=c("sample","variable"),
                           direction="wide")
colnames(all_weights_long)[3:ncol(all_weights_long)] = all_weights$time

# Add pseudocounts and rescale
all_weights_long[,3:ncol(all_weights_long)] = all_weights_long[,3:ncol(all_weights_long)] + 0.001
all_weights_long[,3:ncol(all_weights_long)] = sweep(all_weights_long[,3:ncol(all_weights_long)],2,colSums(all_weights_long[,3:ncol(all_weights_long)], na.rm=TRUE),`/`)

# If there are early and late columns, calculate early and late changes
if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
  all_weights_long$early_late = log2((all_weights_long$`clonal [late]`/(1 - all_weights_long$`clonal [late]`)) / (all_weights_long$`clonal [early]`/(1 - all_weights_long$`clonal [early]`)))
}

# Same for clonal and subclonal
if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
  all_weights_long$clonal_subclonal = log2((all_weights_long$subclonal/(1 - all_weights_long$subclonal)) / (all_weights_long$clonal/(1 - all_weights_long$clonal)))
}


# Bootstrapping for early-late and clonal-subclonal

if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
  changes_1000_el = as.data.frame.matrix(t(replicate(1000, replicateSignatureChanges(multi_list=events, time="el"), simplify="vector")))
}

if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
  changes_1000_cs = as.data.frame.matrix(t(replicate(1000, replicateSignatureChanges(multi_list=events, time="cs"), simplify="vector")))
}
  
# Get confidence intervals for early-late signature changes
all_weights_long = all_weights_long[!is.na(rowSums(all_weights_long[,3:ncol(all_weights_long)])),]

CIs = apply(all_weights_long, 1, function(x) {
  
  sig = x[2] # which signature?
  
  # Early/late confidence intervals
  if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
    el_reps = unlist(changes_1000_el[,sig])
    el_lower = quantile(as.numeric(el_reps), 0.025, na.rm=TRUE)
    el_upper = quantile(as.numeric(el_reps), 0.975, na.rm=TRUE)
  } else { 
    el_lower=NA
    el_upper=NA
  }
 

  # Clonal/subclonal
  if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
    cs_reps = unlist(changes_1000_cs[,sig])
    cs_lower = quantile(as.numeric(cs_reps), 0.025, na.rm=TRUE)
    cs_upper = quantile(as.numeric(cs_reps), 0.975, na.rm=TRUE)  
  } else {
    cs_lower=NA
    cs_upper=NA
  }
  
  c(el_lower, el_upper, cs_lower, cs_upper)
  
})

all_weights_long$el_lCI = CIs[1,]
all_weights_long$el_uCI = CIs[2,]
all_weights_long$cs_lCI = CIs[3,]
all_weights_long$cs_uCI = CIs[4,]

# Add n to all_weights, remove if both 0
n_weights$time = rownames(n_weights)
n_weights_m = melt(n_weights)
n_weights_m = reshape(n_weights_m,
                      timevar="time",
                      idvar=c("sample", "variable"),
                      direction="wide")
n_prop_weights = merge(all_weights_long, n_weights_m)

# Save all results 
write.table(n_prop_weights, paste0(id, "_sig_weights.txt"), sep="\t", quote=FALSE)

# Combine bootstrapping results
if (exists("changes_1000_cs")){
  changes_1000_cs = apply(changes_1000_cs, 2, unlist)
  mode(changes_1000_cs) = "numeric"
  changes_1000_cs = data.frame(changes_1000_cs)
  changes_1000_cs$time = "clonal_subclonal"
} 

if (exists("changes_1000_el")){
  changes_1000_el = apply(changes_1000_el, 2, unlist)
  mode(changes_1000_el) = "numeric"
  changes_1000_el = data.frame(changes_1000_el)
  changes_1000_el$time = "early_late"
}

if (exists("changes_1000_el") & exists("changes_1000_cs")){
  
	print("all")
  all_bootstraps = rbind(changes_1000_cs, changes_1000_el)
  write.table(all_bootstraps, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", quote=FALSE)

} else if(exists("changes_1000_cs") & !exists("changes_1000_el")){
	print("cs")
    write.table(changes_1000_cs, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", quote=FALSE)
  
} else if(exists("changes_1000_el") & !exists("changes_1000_cs")){
	print("el")
    write.table(changes_1000_el, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", quote=FALSE)
  
}


