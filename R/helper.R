#' Weighted Non-negative Sparse Error Recovery
#'
#' Y_\{G x J\} ~= X_\{G x K\} alpha_\{K x J\} + R_\{G x J\}
#'
#' @param Y Gene x Grid (or Spot) cpm-normalized count matrix
#' @param X Gene x Feature matrix: either (1) average expression matrix or (2) NMF feature matrix
#' @param Y.nUMI vector of UMI counts of unnormalized Y
#' @param C tuning parameter
#' @param max.iter maximum number of iterations
#' @param eps convergence tolerance
#'
#' @return
#' \item{loss}{vector of loss function values at each iteration}
#' \item{w}{estimated weight matrix (Gene x Grid)}
#' \item{alpha.hat}{estimated cell type proportion matrix (Feature x Grid)}
#' \item{R.hat}{estimated residual matrix (Gene x Grid)}
#' \item{L.hat}{estimated cellness matrix (Gene x Grid)}
#' @export
wnn.ser <- function(Y, X, Y.nUMI, R.init = NULL, C=2, max.iter = 100, eps = 1e-4, verbose = FALSE){
  
  F2norm <- function(Y) {sqrt(sum(Y^2))}
  thr.op <- function(x,thr) {sign(x)*pmax(abs(x)-thr,0)}
  thresh.l1 <- function(x,thr) {thr.op(x,thr)}
  
  G <- nrow(Y); J <- ncol(Y); K <- ncol(X)
  R <- matrix(0,nrow=G,ncol=J); alpha <- matrix(0,nrow=K,ncol=J)
  mu = prod(dim(Y))/(C*sum(abs(Y))); imu <- 1/mu
  XtX <- crossprod(X)
  
  if(!is.null(R.init)) R <- R.init
  
  i <- 0; obj.cache <- 0; rel.cache <- c()
  
  while(TRUE){
    i <- i+1
    
    # weight update
    if(i==1){
      W = matrix(1, nrow=G, ncol=J)
    } else{
      W = pmax(sqrt(sweep(R, 2, Y.nUMI, "/") * 10^6), 1)
    }
    
    # alpha update
    YR <- Y-R; XtYR <- crossprod(X, YR)
    alpha <- sapply(1:J, function(j) RcppML::nnls(a=XtX, b=as.matrix(XtYR[,j]), fast_nnls=TRUE))
    
    # R update
    wimu <- ( (1/mu) / W )  # (1/mu) * (1/sqrt(Var(R_gj)))
    R <- thresh.l1((Y-X%*%alpha), wimu)
    R[Y==0] <- 0
    R <- pmax(R, 0)
    
    # results
    obj1 <- sum(abs(R/W))
    obj2 <- (mu/2) * F2norm((Y - X%*%alpha - R))^2
    obj <- obj1 + obj2
    
    obj.cache = c(obj.cache, obj)
    rel <- abs(obj - obj.cache[i]) / obj.cache[i]
    rel.cache <- c(rel.cache,rel)
    
    if (verbose) print(c(iter=i, obj=obj))
    if ((i > max.iter) || rel < eps)
      break;
    
  }
  
  result = list(loss=obj.cache[-1], W=W, alpha.hat=alpha, R.hat=R, L.hat=Y-R)
  return(result)
}

#' Naive Non-negative Least Squares (NNLS)
#'
#' Y_\{G x J\} ~= X_\{G x K\} alpha_\{K x J\} + R_\{G x J\}
#'
#' @param Y Gene x Grid (or Spot) cpm-normalized count matrix
#' @param X Gene x Feature matrix: either (1) average expression matrix or (2) NMF feature matrix
#' @param ncores number of cores to use
#' @param chunk.size number of chunk to split
#'
#' @import parallel
#' @import data.table
#'
#' @return
#' \item{alpha.hat}{estimated cell type proportion matrix (Feature x Grid)}
#' \item{resid}{estimated residual value matrix (Gene x Grid)}
#' \item{resid.truncated}{non-negative truncated residual value matrix (Gene x Grid)}
#' 
#' @export
nnls.naive <- function(Y, X, ncores=NULL, chunk.size=NULL){
  
  # Fit NNLS per grids (parallel)
  if(is.null(ncores)) ncores <- detectCores()-1
  if(is.null(chunk.size)) chunk.size = round(ncol(Y) / 100)
  
  chunks = data.table(t(Y))
  chunks[, chunk := floor(.I / chunk.size)]
  chunks = split(chunks, by = "chunk", keep.by = FALSE)
  chunks = lapply(chunks, data.matrix)
  
  cl = makeCluster(ncores)
  clusterExport(cl, c("X"), envir = environment())
  
  f0 <- function(chunk) {
    apply(chunk, 1, function(y) {
      nnls::nnls(A=X, b=y)$x
    })
  }
  
  environment(f0) <- environment(nnls.naive)
  out = parLapply(cl, chunks, f0)
  stopCluster(cl)
  
  
  # results
  alpha <- matrix(unlist(out), nrow=ncol(X), ncol=ncol(Y))
  
  resid <- Y - X %*% alpha
  resid.truncated <- pmax(Y - X %*% alpha, 0)
  
  result = list(alpha.hat=alpha, resid=resid, resid.truncated=resid.truncated)
  return(result)
  
}

#' Modified get_cell_type_info() function from RCTD (Cable et al., Nature Biotechnology, 2022)
#' 
#' @param raw.data a Digital Gene Expression matrix, with gene names as rownames and single cells as columns (barcodes for colnames)
#' @param cell_types a named list of cell type assignment for each cell in \code{raw.data}
#' @param nUMI a named list of total UMI count for each cell in \code{raw.data}
#' @param cell_type_names a named list of cell types
#' 
#' @return
#' \item{cell_type_info}{a list of three elements: (1) \code{cell_type_means} (a data_frame (genes by cell types) for mean normalized expression) (2) \code{cell_type_names} (a list of cell type names) and (3) the number of cell types}
#' 
#' @export
mod.get_cell_type_info <- function(raw.data, cell_types, nUMI, cell_type_names = NULL) {
  if(is.null(cell_type_names))
    cell_type_names = levels(cell_types)
  
  n_cell_types = length(cell_type_names)
  
  get_cell_mean <- function(cell_type) {
    cell_type_data = raw.data[,cell_types == cell_type]
    cell_type_umi = nUMI[cell_types == cell_type]
    normData = sweep(cell_type_data,2,cell_type_umi,`/`)
    return(rowSums(normData) / dim(normData)[2])
  }
  
  cell_type = cell_type_names[1]
  cell_type_means <- data.frame(get_cell_mean(cell_type))
  colnames(cell_type_means)[1] = cell_type
  for (cell_type in cell_type_names[2:length(cell_type_names)]) {
    cell_type_means[cell_type] = get_cell_mean(cell_type)
  }
  return(list(cell_type_means, cell_type_names, n_cell_types))
}

#' Modified get_de_genes() function from RCTD (Cable et al., Nature Biotechnology, 2022)
#' 
#' @param cell_type_info cell type information and profiles of each cell, calculated from the scRNA-seq reference (see \code{\link{get_cell_type_info}})
#' @param Y.count Gene x Grid (or Spot) raw count matrix
#' @param MIN_OBS the minimum number of occurrences of each gene in Y.count (Default: 10)
#' @param fc_thresh minimum \code{log_e} fold change required for a gene (Default: 0.75)
#' @param expr_thresh minimum expression threshold, as normalized expression (Default: 2e-04)
#'
#' @return
#' \item{total_gene_list}{a list of differntially expressed gene names}
#' 
#' @export
mod.get_de_genes <- function(cell_type_info, puck, fc_thresh = 0.75, expr_thresh = 2e-04, MIN_OBS = 10) {
  total_gene_list = c()
  epsilon = 1e-9
  bulk_vec = rowSums(puck)
  gene_list = rownames(cell_type_info[[1]])
  prev_num_genes <- min(length(gene_list), length(names(bulk_vec)))
  if(length(gene_list) == 0)
    stop("get_de_genes: Error: 0 common genes between SpatialRNA and Reference objects. Please check for gene list nonempty intersection.")
  gene_list = gene_list[bulk_vec[gene_list] >= MIN_OBS]
  if(length(gene_list) < 0.1 * prev_num_genes)
    stop("get_de_genes: At least 90% of genes do not match between the SpatialRNA and Reference objects. Please examine this. If this is intended, please remove the missing genes from the Reference object.")
  for(cell_type in cell_type_info[[2]]) {
    if(cell_type_info[[3]] > 2)
      other_mean = rowMeans(cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type])
    else {
      other_mean <- cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type]
      names(other_mean) <- gene_list
    }
    logFC = log(cell_type_info[[1]][gene_list,cell_type] + epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_info[[1]][gene_list,cell_type] > expr_thresh)) #| puck_means[gene_list] > expr_thresh)
    message(paste0("get_de_genes: ", cell_type, " found DE genes: ",length(type_gene_list)))
    total_gene_list = union(total_gene_list, type_gene_list)
  }
  total_gene_list = gene_list[total_gene_list]
  message(paste0("get_de_genes: total DE genes: ",length(total_gene_list)))
  return(total_gene_list)
}

#' Modified remove_platform_effect() function from RCTD (Cable et al., Nature Biotechnology, 2022)
#' 
#' @param Y.count Gene x Grid (or Spot) raw count matrix
#' @param X.count Gene x Cell raw count matrix from reference data
#' @param pseudocount pseudo-counts to be added (Default: 1)
#'
#' @return
#' \item{X.count_corrected}{Gene x Cell platform-effect-normalized count matrix}
#' 
#' @export
mod.remove_platform_effect <- function(
    Y.count,   
    X.count,  
    pseudocount = 1
) {
  # 1) Compute a pseudo-bulk sum for each gene in the ST data
  S_j <- rowSums(Y.count)
  
  # 2) Compute a pseudo-bulk sum for each gene in the reference
  R_j <- rowSums(X.count)
  
  # 3) Estimate a per-gene log-scale factor gamma_j.
  gamma_j <- log(S_j + pseudocount) - log(R_j + pseudocount)
  
  # 4) Rescale the reference by exp(gamma_j) on each gene
  X.count_corrected <- sweep(X.count, 1, exp(gamma_j), `*`)
  
  return(X.count_corrected)
}