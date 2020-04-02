#' @title eval
#'
#' @description function to evaluate predictions from a spectroscopic model
#'
#' @param obs vector of observed values
#' @param pred vector of predicted values
#' @param obj either `cat` for categorical variables of `quant` for quantitative variables
#'
#' @return a set of accuracy measures
#'
#' @export

eval <- function(obs, pred, obj){

  if(obj == 'quant'){

    # mean error
    ME <- round(mean(pred - obs, na.rm = TRUE), digits = 2)

    # root mean square error
    RMSE <-   round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2)

    # Pearson's correlation squared
    r2 <-  round((cor(pred, obs, method = 'spearman', use = 'pairwise.complete.obs')^2), digits = 2)

    # coefficient of determination
    SSE <- sum((pred - obs) ^ 2, na.rm = T)
    SST <- sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
    R2 <- round((1 - SSE/SST), digits = 2)

    # concordance correlation coefficient
    n <- length(pred)
    sdPred <- sd(pred, na.rm = T)
    sdObs <- sd(obs, na.rm = T)
    r <- stats::cor(pred, obs, method = 'pearson', use = 'pairwise.complete.obs')
    # scale shift
    v <- sdPred / sdObs
    sPred2 <- var(pred, na.rm = T) * (n - 1) / n
    sObs2 <- var(obs, na.rm = T) * (n - 1) / n
    # location shift relative to scale
    u <- (mean(pred, na.rm = T) - mean(obs, na.rm = T)) / ((sPred2 * sObs2)^0.25)
    Cb <- ((v + 1 / v + u^2)/2)^-1
    rCb <- r * Cb
    rhoC <- round(rCb, digits = 2)

    # RPD
    sdObs <- sd(obs)
    RMSE_2 <- sqrt(mean((pred - obs)^2))
    RPD <- round((sdObs/RMSE_2), digits = 2)

    # RPIQ
    q25 <- as.numeric(quantile(obs)[2])
    q75 <- as.numeric(quantile(obs)[4])
    iqDist <- q75 - q25
    RMSE <- sqrt(mean((pred - obs)^2))
    rpiq <- round((iqDist/RMSE), digits = 2)

    # return the results
    evalRes <- data.frame(ME = ME, RMSE = RMSE, r2 = r2, R2 = R2, rhoC = rhoC, RPD = RPD, RPIQ = rpiq)
  }
  if (obj =='cat'){

    # overall accuracy
    cm = as.matrix(table(obs = obs, pred = pred))
    n <- length(obs)
    diag = diag(cm)
    OA <- round((sum(diag) / n), digits = 2)

    # Cohens' kappa
    cm = as.matrix(table(obs = obs, pred = pred))
    rowsums = apply(cm, 1, sum)
    colsums = apply(cm, 2, sum)
    n <- length(obs)
    diag = diag(cm)
    accuracy <- sum(diag) / n
    p = rowsums / n
    q = colsums / n
    expAccuracy = sum(p*q)
    kappa = round(((accuracy - expAccuracy) / (1 - expAccuracy)), digits = 2)

    # return the results
    evalRes <- data.frame(OA = OA, kappa = kappa)
  }
  return(evalRes)
}
