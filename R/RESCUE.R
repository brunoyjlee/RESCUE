#' RESCUE
#'
#' Y_\{G x J\} ~= X_\{G x K\} alpha_\{K x J\} + R_\{G x J\}
#'
#' @param Y.count Gene x Grid (or Spot) raw count matrix
#' @param X.count Gene x Cell raw count matrix from reference data
#' @param X.cluster vector of cell type annotations for each cell from reference data
#' @param C.grid vector of tuning parameters; Default: seq(from=1e-6, to=5, length=10)
#' @param seed random seed for generating Poisson surrogate model (Default: 1) 
#' @param max.iter maximum number of iterations (Default: 100)
#' @param eps convergence tolerance (Default: 1e-4)
#' @param platform_norm if true, apply platform effect normalization (Default: TRUE)
#' @param fc_thresh minimum log-fold-change (across cell types) in gene filtering (Default: 0.75) 
#' @param expr_thresh minimum normalized gene expression in gene filtering (Default: 2e-4)
#' @param MIN_OBS minimum total counts per pixel of genes (Default: 10)
#' @param CELL_MIN_INSANCE minimum number of cells requires per cell type (Default: 25)
#' @param ncores number of cores to use; Default: detectCores()-1
#' @param chunk.size number of chunk to split; Default: round(ncol(Y.count) / 100)
#'
#' @import parallel
#' @import data.table
#'
#' @return
#' \item{alpha.hat}{estimated cell type proportion matrix (Feature x Grid)}
#' \item{R.hat}{estimated residual matrix in CPM scale (Gene x Grid)}
#' \item{L.hat}{estimated cellness matrix in CPM scale (Gene x Grid)}
#' \item{R.raw}{estimated residual matrix in raw scale (Gene x Grid)}
#' \item{L.raw}{estimated cellness matrix in raw scale (Gene x Grid)}
#' \item{rel}{vector of relative errors}
#' \item{C.grid}{vector of tuning parameters}
#' 
#' @export

RESCUE <- function(Y.count, X.count, X.cluster, C.grid=NULL, 
                   seed=1, max.iter = 100, eps = 1e-4, platform_norm = TRUE, 
                   fc_thresh = 0.75, expr_thresh = 2e-4, MIN_OBS = 10, CELL_MIN_INSTANCE = 25,
                   ncores=NULL, chunk.size=NULL, verbose=FALSE){
  
  # Pre-processing inputs
  idx.ct <- names(which(table(X.cluster) >= CELL_MIN_INSTANCE))
  if(length(idx.ct)!=length(levels(X.cluster))){
    idx.cell <- which(X.cluster %in% idx.ct)
    X.count <- X.count[,idx.cell]; X.cluster <- factor(as.vector(X.cluster)[idx.cell])  
  }
  X.cluster <- as.factor(X.cluster); X.nUMI <- colSums(X.count)
  
  # Gene filtering
  X.cell_type_info <- mod.get_cell_type_info(X.count, X.cluster, X.nUMI)
  X.get_de_genes <- mod.get_de_genes(cell_type_info=X.cell_type_info, puck=Y.count,
                                    fc_thresh = fc_thresh, MIN_OBS = MIN_OBS, expr_thresh = expr_thresh)
  
  Y.count <- Y.count[X.get_de_genes,]; X.count <- X.count[X.get_de_genes,]
  
  # Check
  if(length(which(colSums(Y.count) == 0))!=0) Y.count <- Y.count[, which(colSums(Y.count)!=0)]
  if(length(which(colSums(X.count) == 0))!=0) X.count <- X.count[, which(colSums(X.count)!=0)]
  Y.nUMI <- colSums(Y.count)
  
  # Platform effect normalization
  if(platform_norm==TRUE){
    X.norm <- mod.remove_platform_effect(Y.count, X.count)  
  } else{
    X.norm <- X.count
  }
  
  # CPM normalization
  X.cpm <- t(t(X.norm) / rowSums(t(X.norm)))*10^6
  X.avg = matrix(0, nrow=nrow(X.norm), ncol=length(levels(X.cluster)))
  for(j in 1:length(levels(X.cluster))){
    X.avg[,j] <- rowMeans(X.cpm[, which(X.cluster==levels(X.cluster)[j])])
  }
  colnames(X.avg) <- levels(X.cluster)
  rownames(X.avg) <- rownames(X.norm)
  
  X <- t((t(X.avg) / rowSums(t(X.avg)))*10^6)
  Y <- t((t(Y.count)/rowSums(t(Y.count)))*10^6)
  
  # Fit NNLS to get alpha.srgt and R.srgt
  X <- as.matrix(X); Y <- as.matrix(Y)
  res.srgt <- nnls.naive(Y=Y, X=X, ncores=ncores, chunk.size=chunk.size)
  gc()
  
  # Construct Poisson surrogate model
  if(verbose==TRUE) print("Constructing Poisson Surrogate Model..")
  alpha.srgt <- res.srgt$alpha.hat
  L.srgt <- X %*% alpha.srgt; R.srgt <- res.srgt$resid.truncated
  Y.srgt <- (L.srgt + R.srgt)
  
  Y.tmp <- t(Y.srgt) / 10^6; rm(alpha.srgt); rm(Y.srgt); gc()
  Y.rescaled <- sweep(Y.tmp, 1, Y.nUMI, "*"); rm(Y.tmp); gc()
  
  set.seed(seed)
  Y.pois.rescaled <- matrix(rpois(prod(dim(Y.rescaled)), Y.rescaled), nrow(Y.rescaled), ncol(Y.rescaled)); rm(Y.rescaled); gc()
  
  # Prevent empty grid if any
  if(length(which(rowSums(Y.pois.rescaled)==0)) !=0){
    idx <- which(rowSums(Y.pois.rescaled)==0)
    Y.pois.rescaled[idx, sample(1:ncol(Y.pois.rescaled), length(idx), replace=TRUE)] <- 1
  }
  
  # Fit WNNSER for tuning parameters over C.grid
  Y.pois.cpm <- t((Y.pois.rescaled / rowSums(Y.pois.rescaled)) * 10^6); rm(Y.pois.rescaled); gc()
  rel.L <- rel.R <- rel <- c()
  if(is.null(C.grid)) C.grid <- seq(from=1e-6, to=5, length=10)
  if(verbose==TRUE) print("Tuning Starts..")  
  for(i in 1:length(C.grid)){
    
    if(verbose==TRUE) print(paste(paste(paste("Tuning Procedure:", i)), "out of", length(C.grid)))
    res.ser <- wnn.ser(Y=Y.pois.cpm, X=X, Y.nUMI=Y.nUMI, R.init=NULL,
                       C = C.grid[i], max.iter=max.iter, eps=eps, verbose=verbose)
    
    alpha.hat <- res.ser$alpha.hat
    R.hat <- res.ser$R.hat
    L.hat <- X%*%alpha.hat
    
    # error metric
    rel.L[i] <- norm(L.srgt - L.hat, type="F") / norm(L.srgt, type="F")  
    rel.R[i] <- norm(R.srgt - R.hat, type="F") / norm(R.srgt, type="F")
    rel[i] <- max(rel.L[i], rel.R[i])
    
    rm(alpha.hat); rm(R.hat); rm(L.hat); gc()
    
  }
  
  # final model
  C.hat <- C.grid[which.min(rel)] 
  if(verbose==TRUE) print(paste("Otimal Tuning Parameter:", C.hat))
  if(verbose==TRUE) print("Fitting the Final Model..")  
  res.ser <- wnn.ser(Y=Y, X=X, Y.nUMI=Y.nUMI, R.init=NULL,
                     C=C.hat, max.iter=max.iter, eps=eps, verbose=verbose)
  
  alpha.hat <- res.ser$alpha.hat
  R.hat <- res.ser$R.hat
  L.hat <- res.ser$L.hat
  
  # Rescale to the raw scale
  R.raw <- (sweep(R.hat/10^6, 2, colSums(Y.count), "*"))
  L.raw <- (sweep(L.hat/10^6, 2, colSums(Y.count), "*"))
  
  if(verbose==TRUE) print("Done!")  
  gc()
  
  result = list(alpha.hat=alpha.hat, Y=Y, X=X,
                R.hat=R.hat, L.hat=L.hat, R.raw=R.raw, L.raw=L.raw, 
                rel=rel, C.grid=C.grid)
  return(result)
}