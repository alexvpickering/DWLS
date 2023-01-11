#' Trim Signature
#'
#' Trim bulk and single-cell data to contain the same genes
#'
#' @param Signature single-cell signature from \code{buildSignatureMatrix}
#' @param bulkData bulk expression matrix (genes x samples)
#'
#' @return list with single-cell and bulk matrices trimmed to common genes
#' @export
#'
trimData <- function(Signature, bulkData) {

  bulkData <- as.matrix(bulkData)
  Genes <- intersect(rownames(Signature), row.names(bulkData))
  B <- bulkData[Genes, ]
  S <- Signature[Genes, ]
  return(list(sig = S, bulk = B))
}


#' Solve OLS
#'
#' solve using OLS, constrained such that cell type numbers>0
#'
#' @export
solveOLS <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  solution <- quadprog::solve.QP(D, d, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution/sum(solution))
}

# return cell number, not proportion do not print output
solveOLSInternal <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  solution <- quadprog::solve.QP(D, d, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

#' Solve Parallel
#'
#' Run solver in parallel
#'
#' @param S single-cell signature matrix
#' @param Bmat matrix of bulk expression signatures
#' @param ncores number of cores
#' @param solver function to use. Default is \code{solveDampenedWLS}.
#'
#' @return results
#' @export
#'
solveParallel <- function(S, Bmat, solver = solveDampenedWLS, ncores = min(ncol(Bmat), parallel::detectCores())) {

  cl <- parallel::makeCluster(ncores)

  parallel::clusterExport(cl=cl, varlist=c("S", "solver"), envir = environment())

  results <- parallel::parApply(cl, Bmat, 2, function(x) solver(S, x))

  parallel::stopCluster(cl)
  return(results)
}

#' Solve Dampened WLS
#'
#' solve using WLS with weights dampened by a certain dampening constant
#'
#' @param S single-cell signature matrix
#' @param B bulk expression profile for one sample
#'
#' @return predicted cell-type proportions
#' @export
#'
solveDampenedWLS <- function(S, B) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- solveOLSInternal(S, B)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- findDampeningConstant(S, B, solution)
  change <- 1
  while (change > 0.01 & iterations < 1000) {
    newsolution <- solveDampenedWLSj(S, B, solution, j)
    # decrease step size for convergence
    solutionAverage <- rowMeans(cbind(newsolution, matrix(solution, nrow = length(solution), ncol = 4)))
    change <- norm(as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  return(solution/sum(solution))
}

# solve WLS given a dampening constant
solveDampenedWLSj <- function(S, B, goldStandard, j) {
  multiplier <- 1 * 2^(j - 1)
  sol <- goldStandard
  ws <- as.vector((1/(S %*% sol))^2)
  wsScaled <- ws/min(ws)
  wsDampened <- wsScaled
  wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened)
  D <- t(S) %*% W %*% S
  d <- t(S) %*% W %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")
  solution <- quadprog::solve.QP(D/sc, d/sc, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# find a dampening constant for the weights using cross-validation
findDampeningConstant <- function(S, B, goldStandard) {
  solutionsSd <- NULL
  # goldStandard is used to define the weights
  sol <- goldStandard
  ws <- as.vector((1/(S %*% sol))^2)
  wsScaled <- ws/min(ws)
  wsScaledMinusInf <- wsScaled
  # ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier) for each, calculate the variance of the dampened weighted solution for a
  # subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier <- 1 * 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)
    for (i in 1:100) {
      set.seed(seeds[i])  #make nondeterministic
      subset <- sample(length(ws), size = length(ws) * 0.5)  #randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit = lm(B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard)/sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }
    solutionsSd <- cbind(solutionsSd, apply(solutions, 1, sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(colMeans(solutionsSd^2))
  return(j)
}

#' @export
solveSVR <- function(S, B) {
  # scaling
  ub = max(c(as.vector(S), B))  #upper bound
  lb = min(c(as.vector(S), B))  #lower bound
  Bs = ((B - lb)/ub) * 2 - 1
  Ss = ((S - lb)/ub) * 2 - 1

  # perform SVR
  model <- e1071::svm(Ss, Bs, nu = 0.5, scale = TRUE, type = "nu-regression", kernel = "linear", cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef < 0)] <- 0
  coef <- as.vector(coef)
  names(coef) <- colnames(S)
  return(coef/sum(coef))
}

#' Get Seurat Markers
#'
#' perform DE analysis using Seurat
#'
#' @param scdata single-cell RNA-seq count matrix
#' @param id cluster assignment
#' @param ... additional arguments to \code{Seurat::FindAllMarkers}
#'
#' @return list of data.frames with markers for each cluster
#' @export
#'
DEAnalysisSeurat <- function(scdata, id, ...) {

  need <- c('Seurat', 'dplyr')
  have <- sapply(need, requireNamespace)
  if (any(!have)) stop('Please install: ', paste0(need[!have], collapse=', '), '.')


  exprObj <- SeuratObject::CreateSeuratObject(counts = as.matrix(scdata), project = "DE")
  SeuratObject::Idents(exprObj) <- as.vector(id)

  all.markers <- Seurat::FindAllMarkers(exprObj, only.pos = TRUE, ...)

  all.markers <- all.markers |>
    dplyr::group_by(.data$cluster) |>
    dplyr::group_split() |>
    as.list()

  res <- list()
  for (df in all.markers) {
    df <- as.data.frame(df)
    row.names(df) <- df$gene
    cluster <- as.character(df$cluster[1])
    res[[cluster]] <- df
  }

  return(res)
}

#' @export
DEAnalysisPresto <- function(scdata, id) {

  need <- c('presto', 'dplyr')
  have <- sapply(need, requireNamespace)
  if (any(!have)) stop('Please install: ', paste0(need[!have], collapse=', '), '.')

  markers <- presto::wilcoxauc(scdata, y = id)

  markers <- markers |>
    dplyr::rename(avg_log2FC = logFC, p_val_adj = padj)  |>
    dplyr::group_by(.data$group) |>
    dplyr::group_split() |>
    as.list()

  res <- list()
  for (df in markers) {
    df <- as.data.frame(df)
    row.names(df) <- df$feature
    group <- as.character(df$group[1])
    res[[group]] <- df
  }

  return(res)
}


#' Build Signature Matrix
#'
#' build signature matrix using genes identified e.g. by DEAnalysis()
#'
#' @param scdata single-cell RNA-seq count matrix
#' @param id cluster assignment
#' @param markers list of data.frames with markers genes for each cluster with columns
#' 'avg_log2FC' and 'p_val_adj' and \code{row.names} as gene symbols.
#' @param diff.cutoff minimum 'avg_log2FC' value
#' @param pval.cutoff minimum 'p_val_adj' value
#'
#' @return signature for to be used for \code{solveDampenedWLS}
#' @export
#'
buildSignatureMatrix <- function(scdata, id, markers, diff.cutoff = 0.5, pval.cutoff = 0.01) {

  numberofGenes <- c()
  for (i in unique(id)) {
    de_group <- markers[[i]]
    DEGenes <- rownames(de_group)[intersect(which(de_group$p_val_adj < pval.cutoff), which(de_group$avg_log2FC > diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.", i, sep = ""), de_group[which(rownames(de_group) %in% DEGenes[nonMir]), ])
    numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
  }

  # need to reduce number of genes for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  # choose matrix with lowest condition number
  conditionNumbers <- c()
  for (G in 50:200) {
    Genes <- c()
    j = 1
    for (i in unique(id)) {
      if (numberofGenes[j] > 0) {
        temp <- paste("cluster_lrTest.table.", i, sep = "")
        temp <- as.name(temp)
        temp <- eval(parse(text = temp))
        temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
        Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
      }
      j = j + 1
    }
    Genes <- unique(Genes)
    # make signature matrix
    ExprSubset <- scdata[Genes, ]
    Sig <- NULL
    for (i in unique(id)) {
      Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
    }
    colnames(Sig) <- unique(id)
    conditionNumbers <- c(conditionNumbers, kappa(Sig))
  }
  G <- which.min(conditionNumbers) + min(49, numberofGenes - 1)  #G is optimal gene number
  #
  Genes <- c()
  j = 1
  for (i in unique(id)) {
    if (numberofGenes[j] > 0) {
      temp <- paste("cluster_lrTest.table.", i, sep = "")
      temp <- as.name(temp)
      temp <- eval(parse(text = temp))
      temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
      Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
    }
    j = j + 1
  }
  Genes <- unique(Genes)
  ExprSubset <- scdata[Genes, ]
  Sig <- NULL
  for (i in unique(id)) {
    Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
  }
  colnames(Sig) <- unique(id)
  return(Sig)
}

