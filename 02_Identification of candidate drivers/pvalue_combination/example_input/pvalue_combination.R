#!Rscript
#Copyright (c) 2015, Institute for Systems Biology
#
#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Author: William Poole
# Ported to R: David L Gibbs
# Email: dgibbs@systemsbiology.org / william.poole@systemsbiology.org / tknijnen@systemsbiology.org
# Created: June 2015

pop.var <- function(x) var(x) * (length(x)-1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))


#Input: raw data vector (of one variable) with no missing samples. May be a list or an array.
#Output Transforemd data vector w.
transformData <- function(data_vector) {
    dvm = mean(data_vector)
    dvsd = pop.sd(data_vector)
    s = (data_vector-dvm)/dvsd
    distr = ecdf(s)
    sapply(s, function(a) -2*log(distr(a)))
}


#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
#       Note: Method does not deal with missing values within the data.
#Output: An m x m matrix of pairwise covariances between transformed raw data vectors
calculateCovariances <- function(data_matrix){
    transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData)
    covar_matrix = cov(transformed_data_matrix)
    covar_matrix
  }


#Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
        df_brown = df_fisher
        c = 1.0
    }
    x = 2.0*sum( -log(p_values) )

    p_brown = pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
    p_fisher = pchisq(df=df_fisher, q=x, lower.tail=FALSE)

    if (extra_info) {
        return(list(P_test=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF=df_brown))
    }
    else {
        return(p_brown)
    }
}

#Input: Either an m x n data matrix with each of m rows representing a variable and each of n columns representing a sample or a covariance table.
#       A vector of m P-values to combine. May be a list or of type numpy.array.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
empiricalBrownsMethod <- function(p_values, data_matrix, covar_matrix, extra_info = FALSE) {
  # inputs must be numeric
    if (missing(covar_matrix)) covar_matrix <- calculateCovariances(data_matrix)
    return(combinePValues(covar_matrix, p_values, extra_info))
}

#

#Input: Either an m x n data matrix with each of m rows representing a variable and each of n columns representing a sample or a covariance table. Should be of type numeric matrix
#       A numeric vector of m P-values to combine.
#Output: A combined P-value using Kost's Method.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
kostsMethod <- function(p_values, data_matrix, covar_matrix, extra_info = FALSE) {
    if (missing(covar_matrix)) covar_matrix <- calculateKostCovariance(data_matrix)
    combinePValues(covar_matrix, p_values, extra_info = extra_info)
}

#Input correlation between two n x n data vectors.
#Output: Kost's approximation of the covariance between the -log cumulative distributions. This is calculated with a cubic polynomial fit.
kostPolyFit <- function(cor) {
    a1 <- 3.263
    a2 <- 0.710
    a3 <- 0.027 #Kost cubic coeficients
    (a1*cor + a2*cor^2 + a3*cor^3)
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numeric matrix.
#       Note: Method does not deal with missing values within the data.
#Output: An m x m matrix of pairwise covariances between the data vectors calculated using Kost's polynomial fit and the pearson correlation function.
calculateKostCovariance <- function(data_matrix) {
    m = nrow(data_matrix)
    covar_matrix = mat.or.vec(m, m)
    for (i in 1:m) {
        for (j in i:m) {
            res0 <- cor.test(data_matrix[i,], data_matrix[j,])
            cor <- res0$estimate
            p_val <- res0$p.value
            covar = kostPolyFit(cor)
            covar_matrix[i, j] = covar
            covar_matrix[j, i] = covar
          }
        }
    return(covar_matrix)
}

#Autors: Grace Tiao, Ziao Lin
#R library for p-value integration

#Remove methods with all NAs in covariance matrices
NA_rm=function(cov_matrix){
	na_matrix=is.na(cov_matrix)
	idx=which(colSums(na_matrix)==dim(cov_matrix)[1])	
	remove=rownames(cov_matrix)[idx]
	cov_matrix=subset(cov_matrix, select = !colnames(cov_matrix) %in% remove)
	cov_matrix=cov_matrix[!rownames(cov_matrix) %in% remove,]
	return(cov_matrix)
}


#Remove interior NAs in a square covariance matrix
NA_rm_interior=function(matrix){
	clean_matrix=na.omit(matrix)
	square_clean_matrix=clean_matrix[,colnames(clean_matrix) %in% rownames(clean_matrix)]
	return(square_clean_matrix)
}

#Remove duplicate rows from dataframe
remove_dups=function(df){
	idx=duplicated(df)
	return(df[!idx,])
}

#QQ plot
qq_plot=function(v, filename, dir){
v = v[order(v)]
n = length(!is.na(v))
qs = 1:n/(n+1)

pdf(file.path(dir, paste(filename, "QQ.pdf", sep= ".")))
plot(-log10(qs), -log10(v), cex=0.5, xlab="Expected p-value (-log10)", ylab="Observed p-value (-log10)", main=filename, col="grey", pch=19)
abline(0,1)
dev.off()
}


empiricalBrownsMethod <- function(covar_matrix, p_values, extra_info = FALSE) {
    return(combinePValues(covar_matrix, p_values, extra_info))
}


calculateCovariances_NA=function(data_matrix){
    transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData_NA)
    covar_matrix = cov(transformed_data_matrix, use="pairwise.complete.obs")
    covar_matrix
  }


empiricalBrownsMethod_NA=function(covar_matrix, p_values, extra_info = FALSE) {
	p_values_clean=p_values[!is.na(p_values)]
	if (length(p_values_clean)<=1) {
		return(rep("NA", times=4))
	} else {
	p_idx=which(p_values %in% p_values_clean)
	covar_matrix_clean=covar_matrix[p_idx,p_idx]
    return(combinePValues(covar_matrix_clean, p_values_clean, extra_info))
    }
}


transformData_NA=function(data_vector) {
	data_vector_clean=data_vector[!is.na(data_vector)]
    dvm = mean(data_vector_clean)
    if (is.na(pop.sd(data_vector_clean))) {
    	return(rep(NA, times=length(data_vector)))
    } else {
    dvsd = pop.sd(data_vector_clean)
    }
    if (dvsd==0) {
    	s = data_vector
    } else {
    s = (data_vector-dvm)/dvsd
    }
    distr = ecdf(s)
    sapply(s, function(a) -2*log(distr(a)))
}


#Function to remove NAs from covariance matrices
#Update function to remove NAs from covariance matrices
outer_edge_NA_rm = function(cov_matrix) {
	idx=which(is.na(cov_matrix)) %% dim(cov_matrix)[1]
	idx[idx==0]=dim(cov_matrix)[1]
	while ((length(idx) > 0) & (1 %in% idx | dim(cov_matrix)[1] %in% idx)) {
	      edge_idx=idx[(idx==1 | idx==dim(cov_matrix)[1])]
	      remove_idx=as.numeric(names(sort(table(edge_idx), decreasing=T)))  # Sort the indices for NAs by frequency, and remove the most frequent one first
	      remove=rownames(cov_matrix)[remove_idx[1]]
	      cov_matrix=subset(cov_matrix, select = !colnames(cov_matrix) %in% remove)
	      cov_matrix=cov_matrix[!rownames(cov_matrix) %in% remove,]
	      idx=which(is.na(cov_matrix)) %% dim(cov_matrix)[1]
	      idx[idx==0]=dim(cov_matrix)[1]
	      }
	return(cov_matrix)
}


BH_adjust=function(pvals){return(p.adjust(pvals, method="BH"))}

#Summarize significant hits per method for automatic method pruning
get_sig=function(pvalues) {
	#Remove duplicate levels
	pvalues=pvalues[!duplicated(pvalues),]
	FDR_vals=apply(pvalues,2,BH_adjust)
	FDR_idx=FDR_vals<.1
	all_sig_genes=colSums(FDR_idx, na.rm=T)
	#df=data.frame(all_sig_genes)
	#names(df)=tissue
	return(all_sig_genes)
}


# Authors: Grace Tiao, Ziao Lin
# Integration of p-values

# Load libraries and arguments
args=commandArgs(TRUE)

# Source the ebm.R package from the following github repo: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM

#source("pvalue_combination/ebm.R")

#source("pvalue_combination/combine_p_values.library.R")
library("reshape")


# Load data
datafile=args[1]
dir=args[2]
observed=read.table(datafile , header=T)

# Remove duplicate rows
observed=remove_dups(observed)

# For each data type, set p-values greater than 1 to NA
observed$ID=as.factor(observed$ID)
observed[observed>1]<-NA

#----------------
# Plot QQs of reported p-values for each method
sort_p=function(pvals){return(as.numeric(pvals[order(pvals)]))}
max_colors=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','black','#b15928', 'lightblue', 'grey', '#8c510a', '#35978f')

df=observed
methods=names(df)[names(df)!="ID"]
df=apply(df[,methods], 2, sort_p)
n=length(df[,1][!is.na(df[,1])])
qs=1:n/(n+1)
obs=df[,1][!is.na(df[,1])]

if (length(methods)==1){
	pdf(file.path(dir, paste("reported_p_values.QQ.pdf", sep=".")))
	plot(-log10(qs), -log10(obs), cex=0.5, xlab="Expected p-value (-log10)", main=paste("Reported P-values", sep=" "), ylab="Observed p-value (-log10)", col="grey", pch=19)
	abline(0,1)
	legend("topleft", pch=19, col="grey", bty="n", legend=colnames(df), cex=.7)
	dev.off()
} else {
	pdf(file.path(dir, paste("reported_p_values.QQ.pdf", sep=".")))
	plot(-log10(qs), -log10(obs), cex=0.5, xlab="Expected p-value (-log10)", main=paste("Reported P-values", sep=" "), ylab="Observed p-value (-log10)", col="grey", pch=19)
	abline(0,1)
	for (i in 2:dim(df)[2]){
		n = length(df[,i][!is.na(df[,i])])
		if (n == 0){
		print (paste('error method', colnames(df)[i], sep = ''))
		}
		qs = 1:n/(n+1)
		obs = df[,i][!is.na(df[,i])]
		points(-log10(qs), -log10(obs), cex=0.5, pch=19, col=max_colors[i])
	}
	legend("topleft", pch=19, bty="n", col=c("grey", max_colors[2:dim(df)[2]]), legend=colnames(df), cex=.7)
	dev.off()
}

#--------------
# Assess covariance

# Cap p-values to 1e-16
observed[observed<1e-16]<-1e-16

# Calculate covariance matrix on observed data (raw p-values are transformed to -2logp)
m=t(data.matrix(observed[,2:dim(observed)[2]]))
observed_cov=calculateCovariances_NA(m)

# Remove methods where all the entries are NA
cov_matrix=observed_cov
cov_matrix=NA_rm(cov_matrix)
#Remove all outer edges with NAs until there are no more NAs on edges; then remove NAs within the matrix with na.omit
n=dim(cov_matrix)[1]
na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
while (na_sum > 0) {
	cov_matrix=outer_edge_NA_rm(cov_matrix)
	n=dim(cov_matrix)[1]
	na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
}
cov_matrix=NA_rm_interior(cov_matrix)
methods=colnames(cov_matrix)

# Combine p-values
observed_tmp=data.matrix(observed[,methods])
observed_brown=data.frame(matrix(unlist(apply(observed_tmp, 1, empiricalBrownsMethod_NA, covar_matrix=cov_matrix, extra_info=T)), ncol=4, byrow=T))
names(observed_brown)=c("Brown_observed", "Fisher_observed", "Brown_Scale_C_observed", "Brown_DF_observed")

evaluate=cbind(observed[,c("ID", methods)], observed_brown)
write.table(evaluate, file.path(dir, paste("combined_p_values.txt", sep=".")), col.names=T, row.names=F, sep="\t", quote=F)


#--------------------
# Summarize significant genes per method (non-integrated)

sig_gene_counts=get_sig(observed[,methods])
median_sig=median(sig_gene_counts)
c=4
upper_thresh=max(median_sig,1)*c

# Remove only inflated methods for now
upper_outliers=names(sig_gene_counts)[sig_gene_counts>upper_thresh]
keep=names(sig_gene_counts)[!names(sig_gene_counts) %in% upper_outliers]


# Prepare report
removal_report=cbind(sig_gene_counts, median_sig, c, upper_thresh, sig_gene_counts>upper_thresh)
colnames(removal_report)[5]=c("removed_for_inflation")
write.table(removal_report, file.path(dir, paste("automatic_method_removal_report.txt", sep=".")), col.names=T, row.names=T, sep="\t", quote=F)

#----------------------
# Automatically remove inflated methods and re-run p-value integration

# Remove methods from observed methods and then subset covariance matrices to match
observed_subset=observed[,keep]
methods=colnames(cov_matrix)[colnames(cov_matrix) %in% names(observed_subset)]

#if (length(methods)<2) {next} #Skip data type if there aren't enough methods to combine

observed_tmp=observed_subset[,methods] # Match (subset) observed data to the methods in the data type
cov_matrix=cov_matrix[methods,methods] # Trim covariance matrix accordingly

# Remove methods where all the entries are NA
cov_matrix=NA_rm(cov_matrix)
n=dim(cov_matrix)[1]
na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
while (na_sum > 0) {
	cov_matrix=outer_edge_NA_rm(cov_matrix)
	n=dim(cov_matrix)[1]
	na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
}
cov_matrix=NA_rm_interior(cov_matrix)
methods=colnames(cov_matrix)[colnames(cov_matrix) %in% names(observed_tmp)]


# Combine p-values
observed_tmp=observed_subset[,methods] #Trim (subset) observed data to match
combined=data.frame(matrix(unlist(apply(observed_tmp, 1, empiricalBrownsMethod_NA, covar_matrix=cov_matrix, extra_info=T)), ncol=4, byrow=T))
names(combined)=c("Brown_observed", "Fisher_observed", "Brown_Scale_C_observed", "Brown_DF_observed")

evaluate=cbind(observed[,c("ID", methods)], combined)
write.table(evaluate, file.path(dir, paste("combined_p_values.automatic_method_removal.txt", sep=".")), col.names=T, row.names=F, sep="\t", quote=F)
