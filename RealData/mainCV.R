####
#Running this file can get the results for Wheat dataset in Section 4.2. 
#'GetR_p' is the main function and you can change N in line 15 to get results for different settings.
#The seed for the Monte-Carlo run is 1-200.
#You need to load the package 'parallel' for distributed computation.
#The results are in the folder named like 'N10.RData'.
####
GetR_p<-function(r){
  library(fdapace)
  Data_x<-as.matrix(read.table('/home/zhh/LR/Realdata/3_22/Wheat/whtspec.txt'))
  Data_x<-scale(Data_x)
  Y<-as.matrix(read.table('/home/zhh/LR/Realdata/3_22/Wheat/protein.txt'))
  Y<-(Y-mean(Y))/sd(Y)
  #Y<-as.matrix(read.table('/home/zhh/LR/Realdata/3_22/Wheat/moist.txt'))
  pts_grid<-seq(0,1,length.out=100);l=100;n=80;N=10;M<-15;out2<-c();out1<-c()
  FPCA_CE<-function (Ly, Lt, optns = list()) {
    ######
    GetCEScores <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {
      
      if (length(lambda) != ncol(phi))
        stop('No of eigenvalues is not the same as the no of eigenfunctions.')
      
      if (is.null(sigma2))
        sigma2 <- 0
      
      Sigma_Y <- fittedCov + diag(sigma2, nrow(phi))
      
      MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)
      ret <- mapply(function(yVec, muphisig){
        rst = GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat, muphisig$Sigma_Yi, verbose=optns$verbose)
        return(rst)
      }, 
      y, MuPhiSig)
      
      return(ret)
    }
    GetMuPhiSig <- function(t, obsGrid, mu, phi, Sigma_Y) {
      #obsGrid <- signif(obsGrid, 14)
      ret <- lapply(t, function(tvec) {
        #ind <- match(signif(tvec, 14), obsGrid)
        if(length(tvec)!=0){
          return(list(muVec=approx(obsGrid,mu,tvec)$y, phiMat=matrix(apply(phi,2,function(phivec){return(approx(obsGrid,phivec,tvec)$y)}),nrow = length(tvec)), 
                      Sigma_Yi=matrix(pracma::interp2(obsGrid,obsGrid,Sigma_Y,as.numeric(as.vector(sapply(tvec,function(x){return(rep(x,length(tvec)))}))),rep(tvec,length(tvec)),method = 'nearest'),length(tvec),length(tvec))))
        }
        else{
          return(list(muVec=numeric(0),phiMat=numeric(0),Sigma_Yi=numeric(0)))
        }
      })
      
      return(ret)
    }
    
    
    GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd=NULL, verbose=FALSE) {
      
      if (length(yVec) == 0) {
        if (verbose) {
          warning('Empty observation found, possibly due to truncation')
        }
        return(list(xiEst=matrix(NA, length(lamVec)), xiVar=matrix(NA, length(lamVec), length(lamVec)), fittedY=matrix(NA, 0, 0)))
      }
      # Do all subscripting stuff in R
      if (!is.null(newyInd)) {    
        if (length(yVec) != 1){ 
          newPhi <- phiMat[newyInd, , drop=FALSE]
          newMu <- muVec[newyInd]
          yVec <- yVec[-newyInd]
          muVec <- muVec[-newyInd]
          phiMat <- phiMat[-newyInd, , drop=FALSE]
          Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop=FALSE]  
          return ( fdapace:::GetIndCEScoresCPPnewInd( yVec, muVec, lamVec, phiMat, Sigma_Yi, newPhi, newMu) )
        } else {   
          # This should be an uncommon scenario
          Lam <- diag(x=lamVec, nrow = length(lamVec))
          LamPhi <- Lam %*% t(phiMat)
          LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
          xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
          xiVar <- Lam - LamPhi %*% t(LamPhiSig) 
          return( list(xiEst=xiEst, xiVar = xiVar, fittedY=NA) )
        }
      } 
      return(fdapace:::GetIndCEScoresCPP( yVec, muVec, lamVec, phiMat, Sigma_Yi) )
      # Unfortunately function overloading is not yet available in Rcpp
      # GetIndCEScoresCPPnewInd and GetIndCEScoresCPP are nearly identical.
      
    }
    ######
    firsttsFPCA <- Sys.time()
    CheckData(Ly, Lt)
    inputData <- fdapace:::HandleNumericsAndNAN(Ly, Lt)
    Ly <- inputData$Ly
    Lt <- inputData$Lt
    optns = SetOptions(Ly, Lt, optns)
    numOfCurves = length(Ly)
    CheckOptions(Lt, optns, numOfCurves)
    if (optns$usergrid == FALSE & optns$useBinnedData != "OFF") {
      BinnedDataset <- fdapace:::GetBinnedDataset(Ly, Lt, optns)
      Ly = BinnedDataset$newy
      Lt = BinnedDataset$newt
      optns[["nRegGrid"]] <- min(optns[["nRegGrid"]], BinnedDataset[["numBins"]])
      inputData$Ly <- Ly
      inputData$Lt <- Lt
    }
    obsGrid = sort(unique(c(unlist(Lt))))
    regGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid)
    outPercent <- optns$outPercent
    buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
    rangeGrid <- range(regGrid)
    minGrid <- rangeGrid[1]
    maxGrid <- rangeGrid[2]
    cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * 
                            outPercent[1] - buff & regGrid < minGrid + diff(rangeGrid) * 
                            outPercent[2] + buff]
    ymat <- fdapace:::List2Mat(Ly, Lt)
    firsttsMu <- Sys.time()
    userMu <- optns$userMu
    if (is.list(userMu) && (length(userMu$mu) == length(userMu$t))) {
      smcObj <- fdapace:::GetUserMeanCurve(optns, obsGrid, regGrid, buff)
      smcObj$muDense = ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
    }
    else if (optns$methodMuCovEst == "smooth") {
      smcObj = GetSmoothedMeanCurve(Ly, Lt, obsGrid, regGrid, 
                                    optns)
    }
    else if (optns$methodMuCovEst == "cross-sectional") {
      smcObj = GetMeanDense(ymat, obsGrid, optns)
    }
    mu <- smcObj$mu
    lasttsMu <- Sys.time()
    firsttsCov <- Sys.time()
    if (!is.null(optns$userCov) && optns$methodMuCovEst != "smooth") {
      scsObj <- GetUserCov(optns, obsGrid, cutRegGrid, buff, 
                           ymat)
    }
    else if (optns$methodMuCovEst == "smooth") {
      scsObj = fdapace:::GetSmoothedCovarSurface(Ly, Lt, mu, obsGrid, 
                                                 regGrid, optns, optns$useBinnedCov)
    }
    else if (optns$methodMuCovEst == "cross-sectional") {
      scsObj = GetCovDense(ymat, mu, optns)
      if (length(obsGrid) != length(cutRegGrid) || !identical(obsGrid, 
                                                              cutRegGrid)) {
        scsObj$smoothCov = ConvertSupport(obsGrid, cutRegGrid, 
                                          Cov = scsObj$smoothCov)
      }
      scsObj$outGrid <- cutRegGrid
    }
    sigma2 <- scsObj[["sigma2"]]
    lasttsCov <- Sys.time()
    firsttsPACE <- Sys.time()
    workGrid <- scsObj$outGrid
    muWork <- ConvertSupport(obsGrid, toGrid = workGrid, mu = smcObj$mu)
    eigObj = fdapace:::GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, 
                                               workGrid, optns, muWork = muWork)
    truncObsGrid <- obsGrid
    if (!all(abs(optns$outPercent - c(0, 1)) < .Machine$double.eps * 
             2)) {
      truncObsGrid <- truncObsGrid[truncObsGrid >= min(workGrid) - 
                                     buff & truncObsGrid <= max(workGrid) + buff]
      tmp <- TruncateObs(Ly, Lt, truncObsGrid)
      Ly <- tmp$Ly
      Lt <- tmp$Lt
    }
    muObs <- ConvertSupport(obsGrid, truncObsGrid, mu = mu)
    phiObs <- ConvertSupport(workGrid, truncObsGrid, phi = eigObj$phi)
    if (optns$methodXi == "CE") {
      CovObs <- ConvertSupport(workGrid, truncObsGrid, Cov = eigObj$fittedCov)
    }
    if (optns$methodXi == "CE") {
      if (optns$methodRho != "vanilla") {
        if (is.null(optns$userRho)) {
          if (length(Ly) > 2048) {
            randIndx <- sample(length(Ly), 2048)
            rho <-fdapace:::GetRho(Ly[randIndx], Lt[randIndx], optns, 
                                   muObs, muWork, truncObsGrid, CovObs, eigObj$lambda, 
                                   phiObs, eigObj$phi, workGrid, sigma2)
          }
          else {
            rho <-fdapace:::GetRho(Ly, Lt, optns, muObs, muWork, 
                                   truncObsGrid, CovObs, eigObj$lambda, phiObs, 
                                   eigObj$phi, workGrid, sigma2)
          }
        }
        else {
          rho = optns$userRho
        }
        sigma2 <- rho
      }
      scoresObj <-GetCEScores(Ly, Lt, optns, muObs, truncObsGrid, 
                              CovObs, eigObj$lambda, phiObs, sigma2)
    }
    else if (optns$methodXi == "IN") {
      scoresObj <- mapply(function(yvec, tvec) GetINScores(yvec, 
                                                           tvec, optns = optns, obsGrid, mu = muObs, lambda = eigObj$lambda, 
                                                           phi = phiObs, sigma2 = sigma2), Ly, Lt)
    }
    if (optns$fitEigenValues) {
      fitLambda <- FitEigenValues(scsObj$rcov, workGrid, eigObj$phi, 
                                  optns$maxK)
    }
    else {
      fitLambda <- NULL
    }
    lasttsPACE <- Sys.time()
    ret <-fdapace:::MakeResultFPCA(optns, smcObj, muObs, scsObj, eigObj, 
                                   inputData = inputData, scoresObj, truncObsGrid, workGrid, 
                                   rho = if (optns$methodRho != "vanilla") 
                                     rho
                                   else 0, fitLambda = fitLambda, timestamps = c(lasttsMu, 
                                                                                 lasttsCov, lasttsPACE, firsttsFPCA, firsttsMu, firsttsCov, 
                                                                                 firsttsPACE))
    if (optns$plot) {
      plot.FPCA(ret)
    }
    return(ret)
  }
  GetSmoothedMeanCurve <- function (y, t, obsGrid, regGrid, optns){
    
    # For the case of binned data one may use a weighted mean response for each time-point.
    # This is not currently implemented. \hat{y}_i = \sum_i w_i y_i where w_i are the
    # same points for common t_is so we have: \hat{y}_i = n_t w_i \bar{y}
    
    userMu = optns$userMu;
    methodBwMu = optns$methodBwMu;
    npoly = 1
    nder = 0 
    userBwMu = optns$userBwMu; 
    kernel = optns$kernel
    
    # If the user provided a mean function use it
    if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
      
      buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
      rangeUser <- range(optns$userMu$t)
      rangeObs <- range(obsGrid)
      if( rangeUser[1] > rangeObs[1] + buff || 
          rangeUser[2] < rangeObs[2] - buff   ) {
        stop('The range defined by the user provided mean does not cover the support of the data.')
      }
      
      mu = spline(userMu$t, userMu$mu, xout= obsGrid)$y;
      muDense = spline(obsGrid,mu, xout=regGrid)$y;
      bw_mu = NULL;
      
      # otherwise if the user provided a mean bandwidth use it to estimate the mean function (below)
    } else {
      if (userBwMu > 0){
        bw_mu = userBwMu;
        #otherwise estimate the mean bandwith via the method selected to estimnate the mean function (below)
      } else {
        if( any(methodBwMu == c('GCV','GMeanAndGCV') )){
          # get the bandwidth using GCV
          bw_mu =  unlist(GCVLwls1D1(yy = y, tt = t, kernel = kernel, npoly = npoly, nder = nder, dataType = optns$dataType) )[1]    
          if ( 0 == length(bw_mu)){ 
            stop('The data is too sparse to estimate a mean function. Get more data!\n')
          }
          # Uncomment to ensure MATLAB compatibility (AdjustBW1 is removed (3-Jun-16); check older versions.)
          # bw_mu = AdjustBW1(kernel=kernel,bopt=bw_mu,npoly=npoly,dataType=optns$dataType,nder=nder)
          # get the geometric mean between the minimum bandwidth and GCV bandwidth to estimnate the mean function (below)         
          if ( methodBwMu == 'GMeanAndGCV') {
            minbw = Minb( unlist(t),2)
            bw_mu = sqrt(minbw*bw_mu);
          } 
        } else {
          # get the bandwidth using CV to estimnate the mean function (below)
          bw_mu = CVLwls1D(y, t, kernel= kernel, npoly=npoly, nder=nder, dataType= optns$dataType, kFolds = optns$kFoldMuCov, 
                           useBW1SE = optns$useBW1SE); 
        }
      }
      # Get the mean function using the bandwith estimated above:
      xin = unlist(t);    
      yin = unlist(y)[order(xin)];
      xin = sort(xin);    
      win = rep(1, length(xin));
      mu = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = obsGrid, win = win)
      muDense = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = regGrid, win = win)
    }  
    
    result <- list( 'mu' = mu, 'muDense'= muDense, 'bw_mu' = bw_mu);
    class(result) <- "SMC"
    return(result)
  }
  CVr<-function(sampX,Y,kfold=5,method='IN',S){
    CV<-rep(0,8)
    CV_In<-rep(0,8)
    CV_split<-rep(0,8)
    CV_pi<-rep(0,8)
    CV_pace<-rep(0,8)
    n<-length(Y)
    N<-length(sampX$Ly[[1]])
    L<-floor(n/kfold)
    X_y<-matrix(0,n,N)
    X_t<-matrix(0,n,N)
    for (i in 1:n) {
      X_y[i,]<-sampX$Ly[[i]]
      X_t[i,]<-sampX$Lt[[i]]
    }
    for (v in 1:kfold) {
      indout<-((v-1)*L+1):(v*L)
      indin<-sort(setdiff(c(1:n),indout))
      X_in<-X_y[indin,];X_out<-X_y[indout,]
      T_in<-X_t[indin,];T_out<-X_t[indout,]
      Y_in<-Y[indin];Y_out<-Y[indout]
      obsLt_in<-list();obsLy_in<-list()
      nin<-length(indin)
      for (i in 1:length(indin)) {
        obsLt_in[[i]]<-T_in[i,]
        obsLy_in[[i]]<-X_in[i,]
      }
      sampX_in<-list(Ly=obsLy_in,Lt=obsLt_in)
      for (m in 1:8) {
        Ly<-sampX_in$Ly;Lt<-sampX_in$Lt
        for (i in 1:nin) {
          Ly[[i]]<-(Y_in[i]-mean(Y_in))*sampX_in$Ly[[i]]
        }
        opt<-SetOptions(Ly,Lt,list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,methodSelectK=m,usergrid=FALSE))
        obsGrid = sort(unique(c(unlist(Lt))))
        g<-GetSmoothedMeanCurve(Ly,Lt,obsGrid,seq(min(obsGrid), max(obsGrid), length.out = opt$nRegGrid),opt)$muDense
        res<-FPCA(sampX_in$Ly,sampX_in$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,methodSelectK=m,usergrid=FALSE))
        Phi_est<-res$phi
        G<-g%*%Phi_est/l
        b<-G*res$lambda^-1
        beta_pi<-Phi_est%*%t(b)
        a_pi<-c()
        for (i in 1:length(Y_in)) {
          tempp<-approx(pts_grid,beta_pi,sampX_in$Lt[[i]],rule = 2)$y
          a_pi[i]<-Y_in[i]-sum(tempp*sampX_in$Ly[[i]])/N
        }
        ##
        Xi_1<-matrix(0,nrow = nin,ncol = m)
        for (i in 1:nin) {
          for (k in 1:m) {
            Xi_1[i,k]<-sum(sampX_in$Ly[[i]]*approx(pts_grid,Phi_est[,k],xout = sampX_in$Lt[[i]],rule = 2)$y)/N
          }
        }
        f1<-lm(Y_in~Xi_1)
        bVec <- as.vector(f1$coef[-1])
        beta_1<-Phi_est%*%bVec
        f_In<-lm(Y_in~res$xiEst)
        beta_In<-Phi_est%*%as.vector(f_In$coef[-1])
        res<-FPCA_CE(sampX_in$Ly,sampX_in$Lt,optns =list(nRegGrid=l,methodXi='CE',maxK=M,FVEthreshold=1,methodSelectK=m,usergrid=FALSE))
        Phi_est<-res$phi
        f_pace<-lm(Y_in~res$xiEst)
        beta_pace<-Phi_est%*%as.vector(f_pace$coef[-1])
        ####
        beta_2s<-matrix(0,nrow = S,ncol = length(beta_1))
        a<-c()
        for (s in 1:S) {
          Index_1<-sort(sample.int(nin,(nin/2)))
          Index_2<-sort(setdiff(c(1:nin),Index_1))
          sampX_1<-list(sampX_in$Lt[Index_1],sampX_in$Ly[Index_1])
          names(sampX_1)<-c('Lt','Ly')
          sampX_2<-list(sampX_in$Lt[Index_2],sampX_in$Ly[Index_2])
          names(sampX_2)<-c('Lt','Ly')
          Y_1<-Y_in[Index_1]
          Y_2<-Y_in[Index_2]
          res1<-FPCA(sampX_1$Ly,sampX_1$Lt,optns =list(nRegGrid=l,methodXi=method,maxK=8,FVEthreshold=1,methodSelectK=m,usergrid=FALSE))
          res2<-FPCA(sampX_2$Ly,sampX_2$Lt,optns =list(nRegGrid=l,methodXi=method,maxK=8,FVEthreshold=1,methodSelectK=m,usergrid=FALSE))
          Phi_est21<-res1$phi
          Phi_est22<-res2$phi
          for (k in 1:m) {
            if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]   }
          }
          Xi_21<-matrix(0,nrow = nin/2,ncol = m)
          Xi_22<-matrix(0,nrow = nin/2,ncol = m)
          for (i in 1:(nin/2)) {
            for (k in 1:m) {
              Xi_21[i,k]<-sum(sampX_in$Ly[[Index_1[i]]]*approx(pts_grid,Phi_est22[,k],xout = sampX_in$Lt[[Index_1[i]]],rule = 2)$y)/N
            }
          }
          for (i in 1:(nin/2)) {
            for (k in 1:m) {
              Xi_22[i,k]<-sum(sampX_in$Ly[[Index_2[i]]]*approx(pts_grid,Phi_est21[,k],xout = sampX_in$Lt[[Index_2[i]]],rule = 2)$y)/N
            }
          }
          Xi_3<-rbind2(Xi_21,Xi_22)
          Y_3<-c(Y_1,Y_2)
          f3<-lm(Y_3~Xi_3)
          a[s]<-f3$coefficients[1]
          bVec3 <- as.vector(f3$coef[-1])
          beta_2s[s,]<-0.5*(Phi_est21+Phi_est22)%*%bVec3
        }
        a_f<-mean(a)
        beta_2<-apply(beta_2s, 2, mean)
        for (i in 1:length(indout)) {
          temp1<-approx(pts_grid,beta_1,T_out[i,],rule = 2)$y
          temp2<-approx(pts_grid,beta_2,T_out[i,],rule = 2)$y
          temp3<-approx(pts_grid,beta_pi,T_out[i,],rule = 2)$y
          temp4<-approx(pts_grid,beta_pace,T_out[i,],rule = 2)$y
          temp5<-approx(pts_grid,beta_In,T_out[i,],rule = 2)$y
          CV[m]<-CV[m]+(Y_out[i]-f1$coefficients[1]-sum(X_out[i,]*temp1)/N)^2
          CV_split[m]<-CV_split[m]+(Y_out[i]-a_f-sum(X_out[i,]*temp2)/N)^2
          CV_pace[m]<-CV_pace[m]+(Y_out[i]-f_pace$coefficients[1]-sum(X_out[i,]*temp4)/N)^2
          CV_In[m]<-CV_In[m]+(Y_out[i]-f_In$coefficients[1]-sum(X_out[i,]*temp5)/N)^2
          CV_pi[m]<-CV_pi[m]+(Y_out[i]-mean(a_pi)-sum(X_out[i,]*temp3)/N)^2
        }
      }
      return(list(CV=which.min(CV),CV_split=which.min(CV_split),CV_In=which.min(CV_In),CV_pi=which.min(CV_pi),CV_pace=which.min(CV_pace)))
    }
  }
  X_train<-Data_x[1:n,]
  Y_train<-Y[1:n]
  X_test<-Data_x[(n+1):100,] 
  Y_test<-Y[(n+1):100] 
  ##
  set.seed(r)
  SampX<-Sparsify(X_train,seq(0,1,length.out=701),N)
  TestX<-Sparsify(X_test,seq(0,1,length.out=701),N)
  m1<-CVr(SampX,Y_train,5,method='IN',S=1)
  m3<-CVr(SampX,Y_train,5,method='IN',S=3)
  m5<-CVr(SampX,Y_train,5,method='IN',S=5)
  ##Estimation
  #plug in
  res<-FPCA(SampX$Ly,SampX$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1, methodSelectK=m1$CV_pi, usergrid=FALSE))
  Phi_est<-res$phi
  Ly<-SampX$Ly;Lt<-SampX$Lt
  for (i in 1:n) {
    Ly[[i]]<-(Y_train[i]-mean(Y_train))*SampX$Ly[[i]]
  }
  opt<-SetOptions(Ly,Lt,list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,methodSelectK=m1$CV_pi, usergrid=FALSE))
  obsGrid = sort(unique(c(unlist(Lt))))
  g<-GetSmoothedMeanCurve(Ly,Lt,obsGrid,seq(min(obsGrid), max(obsGrid), length.out = opt$nRegGrid),opt)$muDense
  G<-g%*%Phi_est/l
  b<-G*res$lambda^-1
  beta_pi<-Phi_est%*%t(b)
  a_pi<-c()
  for (i in 1:length(Y_train)) {
    tempp<-approx(pts_grid,beta_pi,SampX$Lt[[i]],rule = 2)$y
    a_pi[i]<-Y_train[i]-sum(tempp*SampX$Ly[[i]])/N
  }
  ####
  res<-FPCA(SampX$Ly,SampX$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1, methodSelectK=m1$CV_In, usergrid=FALSE))
  Phi_est<-res$phi
  f_In<-lm(Y_train~res$xiEst)
  beta_In<-Phi_est%*%as.vector(f_In$coef[-1])
  ####1
  res<-FPCA(SampX$Ly,SampX$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1, methodSelectK=m1$CV, usergrid=FALSE))
  Phi_est<-res$phi
  Xi_1<-matrix(0,nrow = n,ncol = m1$CV)
  for (i in 1:n) {
    for (k in 1:m1$CV) {
      Xi_1[i,k]<-sum(SampX$Ly[[i]]*approx(pts_grid,Phi_est[,k],xout =SampX$Lt[[i]],rule = 2)$y)/N
    }
  }
  f1<-lm(Y_train~Xi_1)
  bVec <- as.vector(f1$coef[-1])
  beta_1<-Phi_est%*%bVec
  ####pace
  res<-FPCA_CE(SampX$Ly,SampX$Lt,optns =list(nRegGrid=100,methodXi='CE',maxK=M,FVEthreshold=1, methodSelectK=m1$CV_pace, usergrid=FALSE))
  Phi_est<-res$phi
  f_pace<-lm(Y_train~res$xiEst)
  beta_pace<-Phi_est%*%as.vector(f_pace$coef[-1]) 
  #2
  beta_s1<-matrix(0,nrow = 1,ncol = length(beta_1))
  a1<-c()
  for (s in 1:1) {
    Index_1<-sort(sample.int(n,(n/2)))
    Index_2<-sort(setdiff(c(1:n),Index_1))
    SampX_1<-list(SampX$Lt[Index_1],SampX$Ly[Index_1])
    names(SampX_1)<-c('Lt','Ly')
    SampX_2<-list(SampX$Lt[Index_2],SampX$Ly[Index_2])
    names(SampX_2)<-c('Lt','Ly')
    Y_1<-Y_train[Index_1]
    Y_2<-Y_train[Index_2]
    res1<-FPCA(SampX_1$Ly,SampX_1$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m1$CV_split,usergrid=FALSE))
    res2<-FPCA(SampX_2$Ly,SampX_2$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m1$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m1$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m1$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m1$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m1$CV_split) {
        Xi_21[i,k]<-sum( SampX_1$Ly[[i]]*approx(pts_grid,Phi_est22[,k],xout = SampX_1$Lt[[i]],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m1$CV_split) {
        Xi_22[i,k]<-sum(SampX_2$Ly[[i]]*approx(pts_grid,Phi_est21[,k],xout = SampX_2$Lt[[i]],rule = 2)$y)/N
      }
    }
    Xi_3<-rbind2(Xi_21,Xi_22)
    Y_3<-c(Y_1,Y_2)
    f3<-lm(Y_3~Xi_3)
    bVec3 <- as.vector(f3$coef[-1])
    beta_s1[1,]<-0.5*(Phi_est21+Phi_est22)%*%bVec3
    a1<-f3$coefficients[1]
  }
  #3
  beta_s3<-matrix(0,nrow = 3,ncol = length(beta_1))
  a3<-c()
  for (s in 1:3) {
    Index_1<-sort(sample.int(n,(n/2)))
    Index_2<-sort(setdiff(c(1:n),Index_1))
    SampX_1<-list(SampX$Lt[Index_1],SampX$Ly[Index_1])
    names(SampX_1)<-c('Lt','Ly')
    SampX_2<-list(SampX$Lt[Index_2],SampX$Ly[Index_2])
    names(SampX_2)<-c('Lt','Ly')
    Y_1<-Y_train[Index_1]
    Y_2<-Y_train[Index_2]
    res1<-FPCA(SampX_1$Ly,SampX_1$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m3$CV_split,usergrid=FALSE))
    res2<-FPCA(SampX_2$Ly,SampX_2$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m3$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m3$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m3$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m3$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m3$CV_split) {
        Xi_21[i,k]<-sum( SampX_1$Ly[[i]]*approx(pts_grid,Phi_est22[,k],xout = SampX_1$Lt[[i]],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m3$CV_split) {
        Xi_22[i,k]<-sum( SampX_2$Ly[[i]]*approx(pts_grid,Phi_est21[,k],xout = SampX_2$Lt[[i]],rule = 2)$y)/N
      }
    }
    Xi_3<-rbind2(Xi_21,Xi_22)
    Y_3<-c(Y_1,Y_2)
    f3<-lm(Y_3~Xi_3)
    bVec3 <- as.vector(f3$coef[-1])
    beta_s3[s,]<-0.5*(Phi_est21+Phi_est22)%*%bVec3
    a3[s]<-f3$coefficients[1]
  }
  #5
  beta_s5<-matrix(0,nrow = 5,ncol = length(beta_1))
  a5<-c()
  for (s in 1:5) {
    Index_1<-sort(sample.int(n,(n/2)))
    Index_2<-sort(setdiff(c(1:n),Index_1))
    SampX_1<-list(SampX$Lt[Index_1],SampX$Ly[Index_1])
    names(SampX_1)<-c('Lt','Ly')
    SampX_2<-list(SampX$Lt[Index_2],SampX$Ly[Index_2])
    names(SampX_2)<-c('Lt','Ly')
    Y_1<-Y_train[Index_1]
    Y_2<-Y_train[Index_2]
    res1<-FPCA(SampX_1$Ly,SampX_1$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m5$CV_split,usergrid=FALSE))
    res2<-FPCA(SampX_2$Ly,SampX_2$Lt,optns =list(nRegGrid=100,methodXi='IN',maxK=M,FVEthreshold=1,  methodSelectK=m5$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m5$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m5$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m5$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m5$CV_split) {
        Xi_21[i,k]<-sum( SampX_1$Ly[[i]]*approx(pts_grid,Phi_est22[,k],xout = SampX_1$Lt[[i]],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m5$CV_split) {
        Xi_22[i,k]<-sum( SampX_2$Ly[[i]]*approx(pts_grid,Phi_est21[,k],xout = SampX_2$Lt[[i]],rule = 2)$y)/N
      }
    }
    Xi_3<-rbind2(Xi_21,Xi_22)
    Y_3<-c(Y_1,Y_2)
    f3<-lm(Y_3~Xi_3)
    bVec3 <- as.vector(f3$coef[-1])
    beta_s5[s,]<-0.5*(Phi_est21+Phi_est22)%*%bVec3
    a5[s]<-f3$coefficients[1]
  }
  beta_2<-beta_s1[1,]
  beta_3<-apply(beta_s3, 2, mean)
  beta_5<-apply(beta_s5, 2, mean)
  Y_pre1<-c();Y_pre2<-c();Y_pre3<-c();Y_pre5<-c();Y_pi<-c();Y_pace<-c();Y_In<-c()
  for (i in 1:length(Y_test)) {
    temp1<-approx(pts_grid,beta_1,TestX$Lt[[i]],rule = 2)$y
    temp2<-approx(pts_grid,beta_2,TestX$Lt[[i]],rule = 2)$y
    temp3<-approx(pts_grid,beta_3,TestX$Lt[[i]],rule = 2)$y
    tempace<-approx(pts_grid,beta_pace,TestX$Lt[[i]],rule = 2)$y
    temp5<-approx(pts_grid,beta_5,TestX$Lt[[i]],rule = 2)$y
    tempp<-approx(pts_grid,beta_pi,TestX$Lt[[i]],rule = 2)$y
    tempIn<-approx(pts_grid,beta_In,TestX$Lt[[i]],rule = 2)$y
    Y_pre1[i]<-f1$coefficients[1]+sum(temp1*TestX$Ly[[i]])/N
    Y_pre2[i]<-a1+sum(temp2*TestX$Ly[[i]])/N
    Y_pre3[i]<-mean(a3)+sum(temp3*TestX$Ly[[i]])/N
    Y_pre5[i]<-mean(a5)+sum(temp5*TestX$Ly[[i]])/N
    Y_pace[i]<-f_pace$coefficients[1]+sum(tempace*TestX$Ly[[i]])/N
    Y_In[i]<-f_In$coefficients[1]+sum(tempIn*TestX$Ly[[i]])/N
    Y_pi[i]<-mean(a_pi)+sum(tempp*TestX$Ly[[i]])/N
  }
  out2[1]<-sum((Y_test-Y_pre1)^2)/length(Y_test)
  out2[2]<-sum((Y_test -Y_pre2)^2)/length(Y_test)
  out2[3]<-sum((Y_test -Y_pre3)^2)/length(Y_test)
  out2[4]<-sum((Y_test -Y_pre5)^2)/length(Y_test)
  out2[5]<-sum((Y_test -Y_pi)^2)/length(Y_test)
  out2[6]<-sum((Y_test -Y_pace)^2)/length(Y_test)
  out2[7]<-sum((Y_test -Y_In)^2)/length(Y_test)
  out1[1]<-sum(abs(Y_test-Y_pre1))/length(Y_test)
  out1[2]<-sum(abs(Y_test -Y_pre2))/length(Y_test)
  out1[3]<-sum(abs(Y_test -Y_pre3))/length(Y_test)
  out1[4]<-sum(abs(Y_test -Y_pre5))/length(Y_test)
  out1[5]<-sum(abs(Y_test -Y_pi))/length(Y_test)
  out1[6]<-sum(abs(Y_test -Y_pace))/length(Y_test)
  out1[7]<-sum(abs(Y_test -Y_In))/length(Y_test)
  return(list(out1=out1,out2=out2))
}
####
cl <- makeCluster(getOption("cl.cores", 100));
system.time({
  res_p<-parLapply(cl, 1:200,GetR_p)
  stopCluster(cl)
})
####
Outp<-matrix(0,nrow = 200,ncol = 7)
for (i in 1:200) {
  Outp[i,]<-res_p[[i]]$out2
}
apply(Outp, 2, mean)
apply(Outp, 2, sd)
