#' @title cssfunction
#'
#' @description Assesing the adequate calibration set size for infrared spectroscopy
#'
#' @param S A matrix of the socres of the principal components
#'
#' @return css
#'
#' @export


#------------------------------- Info ------------------------------------------
# Description: Function for assesing the adequate calibration set size for
#              - Kennard-Stone Sampling
#              - K-means Sampling
#              - conditioned Latin hypercube Sampling
#
# Authors:     Leo Ramirez-Lopez & Alexandre Wadoux
#              ramirez-lopez.l@buchi.com; alexandre.wadoux@wur.nl
#
# Date:        Jun 2017
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Arguments
#   S:           A matrix of the socres of the principal components
#   k:           A vector containing the sample set sizes to be evaluated
#   method:      The sampleing algorithm. Options are:
#                - "kss" (Kennard-Stone Sampling)
#                - "kms" (K-means Sampling). The default
#                - "clhs" (conditioned Latin hypercube Sampling)
#   repetitions: The number of times that the sampling must be carried out for
#                for each sample size to be evaluated. The results of the
#                final msd is the average of the ones obtained at each iteration.
#                Note that since the "kss" method is deterministic and always return the
#                same results, there is no need for repetitions.
#   n:           The number of equally spaced points at which the probability densities
#                are to be estimated (see density function of the package stats).
#   from, to:    A vector of the left and right-most points of the grid at which the
#                densities are to be estimated. Default is the minimums and maximums
#                of the variables in S.
#   bw:          A vector containing the the smoothing bandwidth to be use for the
#                probability densities (see density function of the package stats).
#   ...:         arguments to be passed to the calibration sampling algorithms, i.e.
#                additional aruments to be used for the clhs, kenStone or naes functions
#                which run inside this function

# Returned object:
# A table with the following columns:
# - css: the sample set size (k)
# - msd
# - msd_sd: the standard deviation of the msd for all the repetitions (does not apply to "kss" since it
# always return the same results)

# Details:
#  This function works by comparing the probability density function (pdf) of the population and the pdf
#  of the sample set in order to asses the representativeness of the sample set.
#  See Ramirez-Lopez,  et al. (2014) for more details.

# References:
#  Ramirez-Lopez, L., Schmidt, K., Behrens, T., van Wesemael, B., Dematt?, J. A., Scholten, T. (2014).
#  Sampling optimal calibration sets in soil infrared spectroscopy. Geoderma, 226, 140-150.

css <- function(S, k, method = "kms", repetitions = 10, n = 512, from, to, bw, ...){
  requireNamespace('clhs')
  requireNamespace('matrixStats')

  if(missing(from)){
    min.sc <- matrixStats::colMins(S)
  }else{
    min.sc <- from
  }
  if(missing(from)){
    max.sc <- matrixStats::colMaxs(S)
  }else{
    max.sc <- to
  }

  if(length(min.sc) != ncol(S))
    stop("Argument 'from' needs to be a vector with the same number of variables (columns) as in S ")

  if(length(max.sc) != ncol(S))
    stop("Argument 'to' needs to be a vector with the same number of variables (columns) as in S ")

  if(missing(bw)){
    d.bandwidths <- apply(S, 2, bw.nrd0)
  }else{
    d.bandwidths <- bw
  }

  if(length(d.bandwidths) != ncol(S))
    stop("Argument 'bw' needs to be a vector with the same number of variables (columns) as in S ")


  ## matrix where the density values will be stored
  sc.dens <- NULL
  for(i in 1:length(min.sc)){
    i.sc.dens <- data.frame(x = seq(min.sc[i], max.sc[i], length = n),
                            densc = rep(NA, n), pc = paste("PC-", i, sep = ""))
    sc.dens <- rbind(sc.dens, i.sc.dens)
  }

  ## estimate the density distribution of each variable
  names(d.bandwidths) <- colnames(S)
  for(i in 1:length(min.sc)){
    idsty <- density(S[,i],
                     bw = "nrd0",
                     n = n, from = min.sc[i], to = max.sc[i],
                     kernel = "gaussian")
    sc.dens[sc.dens$pc == paste("PC-", i, sep = ""),"densc"] <- idsty$y
    #d.bandwidths[i] <- idsty$bw
  }

  ## --- 3. Define the different sample set sizes  ----

  ## --- 4. Sample with the specified algorithm ----
  ## for each sample size the sampling is repeated 10 times and the
  ## differences between the density distribution of the whole set is
  ## compared against the density distribution of the sample set
  ##
  results.ss <-  data.frame(css = k,
                            msd = rep(NA, length(k))
                            #mndiff = rep(NA, length(k)),
                            #sddiff = rep(NA, length(k))
  )
  if(method == "kss" & repetitions > 1){
    warning("For Kenard-Stone Sampling repetitions are not necessary. Only one repetition will be executed")
    repetitions <- 1
  }

  for(i in 1:repetitions){
    results.ss[,-1] <- NA
    fn <- paste(method, "_temp_results_rep", i,".txt", sep = "")

    # if(fn %in% list.files())
    # {
    #   results.kms <- read.table(fn, header = T, sep = "\t")
    # }

    #iter.p <- 1 + sum(rowSums(!is.na(results.kms)) == ncol(results.kms))

    for(j in 1:length(k)){
      set.seed(j)

      if(method == "kms"){
        i.calidx <- prospectr::naes(X = S,
                         k = k[j],
                         method = 0,
                         .center = FALSE,
                         .scale = FALSE, ...)$model
      }

      if(method == "kss"){
        i.calidx <- prospectr::kenStone(X = S,
                             k = k[j],
                             metric = "mahal",
                             .center = FALSE,
                             .scale = FALSE, ...)$model
      }

      if(method == "clhs"){
        i.calidx <- clhs::clhs(as.data.frame(S),
                         size = k[j],
                         simple = TRUE,
                         progress = FALSE, ...)
      }



      m.sc.dens <- msd.sc  <- sc.dens
      for(m in 1:length(min.sc)){

        ## use the same bandwidth (bw) as in the whole set of candidates
        slc <- sc.dens$pc == paste("PC-", m, sep = "")
        m.sc.dens[slc,"densc"] <- density(S[i.calidx, m],
                                          bw = d.bandwidths[m],
                                          n = n, from = min.sc[m], to = max.sc[m],
                                          kernel = "gaussian")$y
      }
      results.ss$msd[j] <- mean((m.sc.dens$densc - sc.dens$densc)^2, na.rm = T)
      #results.kms$mndiff[i] <- mean(abs(colMeans(pcaall$scores.std[i.calidx,])))
      #results.kms$sddiff[i] <- mean(abs(colSds(pcaall$scores.std[i.calidx,]) - 1))

      ## write the results to a table
      if(method == "kss"){
        write.table(results.ss,
                    file = paste(method, "_final_results_.txt", sep = ""),
                    row.names = FALSE, sep = "\t")
        final.ss <- results.ss

      }else{
        write.table(results.ss,
                    file = fn,
                    row.names = FALSE, sep = "\t")
      }
      #print(results.ss[1:i,])
    }
  }

  if(method != "kss"){
    ## --- 5. Read the iteration results from the generated files and compute the mean of the iterations ----
    nmsreps <- paste(method, "_temp_results_rep", 1:repetitions, ".txt", sep = "")
    final.ss <- 0
    for(i in nmsreps){
      iter <- which(i == nmsreps)
      results.ss <- read.table(i, header = T, sep = "\t")
      #results.kms$mndiff <- abs(results.ss$mndiff)
      final.ss <- final.ss + results.ss
      if(iter == length(nmsreps))
      {
        final.ss <- final.ss/iter
      }
    }

    ## --- 6. Read the iteration results from the generated files and compute the standard deviation of the iterations ----
    final.ss_sd <- 0
    for(i in nmsreps){
      iter <- which(i == nmsreps)
      results.ss <- read.table(i, header = T, sep = "\t")
      final.ss_sd <- (results.ss - final.ss_sd)^2
      if(iter == length(nmsreps))
      {
        final.ss_sd <- (final.ss_sd/iter)^0.5
      }
    }
    if(repetitions > 2)
      final.ss  <- data.frame(final.ss, msd_sd = final.ss_sd[,2])
    write.table(final.ss, file = paste(method, "_final_results_.txt", sep = ""), sep = "\t", row.names = FALSE)
  }
  return(final.ss)
}

