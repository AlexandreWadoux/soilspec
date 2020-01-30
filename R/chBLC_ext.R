#' @title chBLC_ext
#'
#' @description Fitting a Convex hull function for regions of interest
#'
#' @param spectra A matrix of spectra
#'
#' @return the fitted convex hull
#'
#' @export

#Convex hull function used for regions of interest work such as detecting secondary clay minerals from vis-NIR spectra
chBLC_ext<- function(spectra, lower, upper, ...){
  interval<- c(1:ncol(spectra))
  hull_spectra<- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
  tempSpect= as.matrix(spectra)
  data1 <- sortedXyData(interval, tempSpect)

  ## calculate convex hull
  c_hull <- chull(data1)
  ## get the appropriate region
  cc<-c_hull
  xs<-which(cc == 1)
  ccd<- cc[which(cc == 1):length(cc)]
  ccd.m<- matrix(NA, nrow=1, ncol=length(ccd))
  for(f in 1:(length(ccd)-1)){
    if(ccd[f] < ccd[f+1]){ccd.m[1,f]<- ccd[f]}
    else {ccd.m[1,f]<- ccd[f]
    break}}
  if(f == (length(ccd)-1)){ccd.m[1,ncol(ccd.m)]<- ccd[length(ccd)]}
  tempX<- c(ccd.m)[!is.na(c(ccd.m))]
  #Left side
  if(c_hull[1] > c_hull[length(c_hull)]){lss<- cc[1:(xs-1)]
  slss<- lss[which(lss >= lss[1])]
  c_hull<- sort(c(tempX,slss))} else {c_hull<-tempX}

  ## calculate linear approximation between hull points
  linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
  ## calculate the deviation from the convex hull
  hull_spectra[1,] <- 1- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])
  colnames(hull_spectra) <- colnames(spectra)
  # construct polygon
  aaa<- hull_spectra
  aaa[,c(1,ncol(aaa))]<- 1.0001
  bbb<- cbind(seq(lower,upper,by =1),t(aaa))
  ccc<- cbind(seq(lower,upper,by =1),rep(c(1.0001),each=length(seq(lower,upper,by=1))))
  ddd<- as.matrix(rbind(bbb,ccc))
  retval<- list(wave=seq(lower,upper,by=1),c.hull=hull_spectra, raw.spec=data1$y,continuum=linear_approx[[2]] ,polygon=ddd)
  return(retval)}
