#' 功能分析gene select 与 DE 的重合个数
#' param:data  矩阵  row:gene    col:cell
#' pbmc.markers 差异表达基因, seurat格式
#' num  挑选的前num个gene
DECompare <- function(data,pbmc.markers,num){
  #dim(data)
  message(paste("first   data: cell ",dim(data)[2]," gene ", dim(data)[1]))
  #message(paste("pbmc.markers: cell ",dim(pbmc.markers)[2]," gene ", pbmc.markers$gene))
  selectgene <- rownames(data)
  if(dim(data)[1] < num)
  {
    message(paste('当前gene个数小于指定值: ', dim(data)[1]))
    DEnumber <- intersect(selectgene,pbmc.markers$gene)
  }
  else{
    DEnumber <- intersect(selectgene[1:num],pbmc.markers$gene)
  }
  return(length(DEnumber))
}

#' Seurat中的聚类方法
#' return 返回标签
Seurat_cluster <- function(data){
  library("Seurat")
  rownames(data) <- paste0("GENE",c(1:dim(data)[1]))
  colnames(data) <- paste0("CELL",c(1:dim(data)[2]))
  seurat_obj <- CreateSeuratObject(data, project = "SEURAT")
  all.genes <- rownames(data)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj,features = rownames(data))
  seurat_obj <- FindNeighbors(seurat_obj,dims = 1:10)
  seurat_obj <- FindClusters(object = seurat_obj)
  return(as.numeric(Idents(seurat_obj)))
}


#非负核自动编码器
scDHA_FS <- function(data = data, k = NULL, K = 3, n = 5000, ncores = 15L, gen_fil = T,sample.prob = NULL, seed = NULL) {
  set.seed(seed)

  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-6
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)

  lr <- c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > n))
  n <- ifelse(gen_fil, min(n, ncol(data)), ncol(data))
  #print(n)
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3

  intermediate_dim = 64
  latent_dim = 15
  if (gen_fil) {
    data.list <- lapply(seq(3), function(i) {
      if(!is.null(seed)) set.seed((seed+i))
      if(nrow(data) > 2000)
      {
        ind <- sample.int(nrow(data), 2000, replace = F, prob = sample.prob)
      } else {
        ind <- seq(nrow(data))
      }
      count <- counts(data)
      #print(typeof(data))
      data.tmp <- as.matrix(count[ind,])
      data.tmp
    })
    print(typeof(data.list))
    or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)
    #or <- data.list

    keep <- which(or > quantile(or,(1-min(n,original_dim)/original_dim)))

    da <- data[,keep]
    original_dim_reduce <- ncol(da)
    or.da <- da
  } else {
    da <- data
    original_dim_reduce <- ncol(da)
    or.da <- da
  }
  return(or.da)
}

gene.filtering <- function(data.list, original_dim, batch_size, ncores.ind, wdecay, seed)
{
  or <- list()
  cl <- parallel::makeCluster(3, outfile = "/dev/null")
  registerDoParallel(cl, cores = 3)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
    library(tensorflow)
  })
  or <- foreach(i = seq(3)) %dopar%
    {
      if (is.null(seed))
      {
        config <- list()
        config$intra_op_parallelism_threads <- ncores.ind
        config$inter_op_parallelism_threads <- ncores.ind
        session_conf <- do.call(tf$ConfigProto, config)
        sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
        keras::k_set_session(session = sess)
      } else {
        set.seed((seed+i))
        use_session_with_seed((seed+i))
      }

      data.tmp <- data.list[[i]]
      batch_size <-round(nrow(data.tmp)/50)

      x <- keras::layer_input(shape = c(original_dim))
      h <- keras::layer_dense(x, 50, kernel_constraint = keras::constraint_nonneg())
      x_decoded_mean <- keras::layer_dense(h, original_dim)
      vae <- keras::keras_model(x, x_decoded_mean)
      magrittr::`%>%`(vae,keras::compile(optimizer =   tf$contrib$opt$AdamWOptimizer(wdecay,1e-3), loss = 'mse'))


      his <- magrittr::`%>%`(vae,keras::fit(
        data.tmp, data.tmp,
        shuffle = TRUE,
        epochs = 10,
        batch_size = batch_size,
        verbose = 0
      ))

      W <- keras::get_weights(keras::get_layer(vae, index = 2))[[1]]
      Wsd <- matrixStats::rowSds(W)
      Wsd[is.na(Wsd)] <- 0
      Wsd <- (Wsd-min(Wsd))/(max(Wsd)-min(Wsd))
      Wsd
    }

  parallel::stopCluster(cl)
  or <- matrixStats::rowMeans2(as.matrix(data.frame(or)))
  or
}


cal_scDHA_FS <- function(data){
  library(SingleCellExperiment)
  message("数据输入：行(gene), 列(细胞)")
  message(paste("检查：row:",dim(data)[1], ", col:", dim(data)[2]))
  counts_matrix <- as.matrix(data)
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
  scdha <- scDHA_FS(sce)
  count <- assay(sce, "counts")
  count
}

Zeroradio <- function(DATA){
  data <- as.matrix(DATA)
  czero <- apply(data, 1, function(x){sum(x == 0)})
  zeronum <- sum(czero)
  r <- zeronum/(length(data))
  r
}
