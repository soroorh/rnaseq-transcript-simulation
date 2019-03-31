
# helper functions -------
txAbundanceByReplicatePerGroup <- function(ntx, nReplicate, shape, scale) matrix(rgamma(ntx*nReplicate, shape = shape, scale = scale), ntx, nReplicate)
as.list.matrix <- function(x) lapply(seq_len(ncol(x)), function(i) x[,i])



# main function -----
## TO DO:
## re-write for a DGEList object?
## argument compatibility checks
## documentation
## author contribution

streamSimReads <- function(aveExpr,dispersions, txEffectiveLength,
                           fc.mat, nReplicates,  
                           nCores = 2, seed,
                           transcript.file, out.sample.size, output.prefix,
                           simplify.transcript.names,
                           read.length = 75,
                           ... ){
  # aveExpr: a named vector of average expressions - deprecated
  # dispersions: named vector of dispersions - deprectaed
  # fc.mat: nTranscript x nGroup matrix of fold-changes of DE transcripts
  # nReplicates: a vector the same length as ncol(fc.mat) specifying the number of biological replicates per group
  # out.sample.size: a vector same length as ncol(fc.mat) of library sizes
 
  
  
  g <- ncol(fc.mat)
  ntx <- nrow(fc.mat)
  
  mu0 <- matrix(rep(aveExpr,g), ncol=g, byrow = FALSE)
  Disp <- matrix(rep(dispersions,g), ncol=g, byrow = FALSE)
  
  mu0 <- mu0 * fc.mat
  
  shape <- 1/Disp
  scale <- mu0/shape
  
  print("Simulating TPM values")
  set.seed(seed = seed)
  
  
  # returns a list of matrices. Each matrix is a transcript x replicate matrix. number of elements of the list == number of groups
  gamma_abun <- mapply(txAbundanceByReplicatePerGroup, rep(ntx, g),
                       nReplicates, as.list.matrix(shape),
                       as.list.matrix(scale),
                       SIMPLIFY = FALSE)

  # transcript x sample matrix of abundance values  sampled from a gamma distribution
  mu <- do.call(cbind, gamma_abun)


  tpm <- mu * read.length * 1e6 / as.numeric(txEffectiveLength)
  tpm <- tpm / matrix(colSums(tpm)/1e6, dim(tpm)[1], dim(tpm)[2], byrow = TRUE)
  contigFile <- Rsubread::scanFasta(transcript.file, simplify.transcript.names = simplify.transcript.names)
  rownames(tpm) <- contigFile$TranscriptID
  
  
  
  print("Parallelizing read simulation")
  # start of parallelization using foreach


  # nCores <- max(nCores,sum(nReplicates))
  # print(paste(nCores, "cores detected"))
  # doParallel::registerDoParallel(nCores)

  
  do.i <- function(i, output.prefix, out.sample.size, 
                   tpm, transcript.file,
                   ...){

    count <- Rsubread::simReads(expression.levels = data.frame(
      TranscriptID=as.character(rownames(tpm)), 
      TPM= tpm[,i]),
                      transcript.file = transcript.file, 
                       output.prefix = output.prefix,
                       out.sample.size = as.integer(out.sample.size),
                      ...)
  
  }
  Nsamples <- ncol(tpm)
  parallel::mcmapply(do.i, 1:Nsamples, 
                             paste(output.prefix,"sample", 1:Nsamples, sep=""),
                             out.sample.size= out.sample.size,
           MoreArgs = list(transcript.file= transcript.file, tpm = tpm,
                           ...),
           mc.cores = nCores,
           mc.preschedule = FALSE
           )

  print("Simulations completed")
  tpm
}


