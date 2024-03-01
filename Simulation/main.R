####
#Running this file can get the results in Section 4.1. 
#'GetResult' is the main function and you can change the n and N in line 11 to get results for different settings.
#The seed for the Monte-Carlo run is 1-200.
#You need to load the package 'parallel' for distributed computation.
#The results are in the folder named like 'N10_n100.RData'.
####
GetResult<-function(r){
  set.seed(r)
  library(fdapace)
  n<-100;N=10;K<-50;M<-8;S<-5;l<-100
  ####
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
            rho <- GetRho(Ly[randIndx], Lt[randIndx], optns, 
                          muObs, muWork, truncObsGrid, CovObs, eigObj$lambda, 
                          phiObs, eigObj$phi, workGrid, sigma2)
          }
          else {
            rho <- GetRho(Ly, Lt, optns, muObs, muWork, 
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
    L<-n/kfold
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
        opt<-SetOptions(Ly,Lt,list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,methodSelectK=m,userSigma2=0.01,usergrid=FALSE))
        obsGrid = sort(unique(c(unlist(Lt))))
        g<-GetSmoothedMeanCurve(Ly,Lt,obsGrid,seq(min(obsGrid), max(obsGrid), length.out = opt$nRegGrid),opt)$muDense
        res<-FPCA(sampX_in$Ly,sampX_in$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m,userSigma2=0.01,usergrid=FALSE))
        Phi_est<-res$phi
        G<-g%*%Phi_est/l
        b<-G*res$lambda^-1
        beta_pi<-Phi_est%*%t(b)
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
        res<-FPCA_CE(sampX_in$Ly,sampX_in$Lt,optns =list(nRegGrid=l,methodXi='CE',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m,userSigma2=0.01,usergrid=FALSE))
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
          res1<-FPCA(sampX_1$Ly,sampX_1$Lt,optns =list(nRegGrid=l,methodXi=method,maxK=8,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m,usergrid=FALSE))
          res2<-FPCA(sampX_2$Ly,sampX_2$Lt,optns =list(nRegGrid=l,methodXi=method,maxK=8,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m,usergrid=FALSE))
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
          CV_pi[m]<-CV_pi[m]+(Y_out[i]-sum(X_out[i,]*temp3)/N)^2
        }
      }
      return(list(CV=which.min(CV),CV_split=which.min(CV_split),CV_In=which.min(CV_In),CV_pi=which.min(CV_pi),CV_pace=which.min(CV_pace)))
    }
  }
  ####generate traning data
  out1<-c()
  out2<-c()
  pts_obs<-matrix(0,nrow = n,ncol = N)
  pts_grid<-seq(0,1,length.out = l)
  basis_grid<-CreateBasis(K,pts_grid,type = 'cos') 
  for (i in 1:n) {
    pts_obs[i,]<-sort(runif(N,0,1))
  }
  Mu<-list(pts_grid,rep(0,l) )
  names(Mu)<-c('t','mu')
  X_grid<-matrix(0, nrow = n, ncol = N)
  Xi_true<-matrix(0, nrow = n, ncol = K)
  b<-matrix(0, nrow = K, ncol = 1)
  for (k in 1:K) {
    if(k==1){b[k]=1}
    else{b[k]<-4*k^{-3}}
  }
  beta_true<-basis_grid%*%b
  #plot(beta_true,type = 'l')
  for (i in 1:n) {
    for (k in 1:K) {
      Xi_true[i,k]<-(1/k)*rnorm(1,0,1)
    }
  }
  X_true=Xi_true%*%t(basis_grid)
  #X_grid=X_true+matrix(0.1*rnorm(n*l), nrow = n, ncol = l)
  for (i in 1:n) {
    phi_i<-CreateBasis(K,pts_obs[i,],type = 'cos') 
    X_grid[i,]<-phi_i%*%Xi_true[i,]+0.1*rnorm(N,0,1)
  }
  obsLt<-list();obsLy<-list()
  for (i in 1:n) {
    obsLt[[i]]<-pts_obs[i,];
    obsLy[[i]]<-X_grid[i,]
  }
  sampX<-list(Ly=obsLy,Lt=obsLt)
  Y<-c()
  for (i in 1:n) {
    Y[i]<-trapzRcpp(pts_grid,beta_true*X_true[i,])+0.5*rnorm(1,0,1)
  }
  ####generate test data
  X_test<-matrix(0, nrow = 2*n, ncol = N)
  X_test_grid<-matrix(0, nrow = 2*n, ncol = l)
  test_obs<-matrix(0,nrow = 2*n,ncol = N)
  for (i in 1:(2*n)) {
    test_obs[i,]<-sort(runif(N,0,1))
  }
  Xi_test<-matrix(0,nrow = 2*n,ncol = K)
  for (i in 1:(2*n)) {
    for (k in 1:K) {
      Xi_test[i,k]<-(1/k)*rnorm(1,0,1)
    }
  }
  for (i in 1:(2*n)) {
    phi_i<-CreateBasis(K,test_obs[i,],type = 'cos') 
    X_test[i,]<-phi_i%*%Xi_test[i,]+0.1*rnorm(N,0,1)
  }
  X_test_grid=Xi_test%*%t(basis_grid)
  testLt<-list();testLy<-list()
  for (i in 1:(2*n)) {
    testLt[[i]]<-test_obs[i,];
    testLy[[i]]<-X_test[i,]
  }
  testX<-list(Ly=testLy,Lt=testLt)
  Y_test<-c()
  for (i in 1:(2*n) ){
    Y_test[i]<-trapzRcpp(pts_grid,beta_true*X_test_grid[i,])+0.5*rnorm(1,0,1)
  }
  ####
  m1<-CVr(sampX,Y,5,method='IN',S=1)
  m3<-CVr(sampX,Y,5,method='IN',S=3)
  m5<-CVr(sampX,Y,5,method='IN',S=5)
  ####Estimation 
  
  #plug in
  res<-FPCA(sampX$Ly,sampX$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m1$CV_pi,userSigma2=0.01,usergrid=FALSE))
  Phi_est<-res$phi
  Ly<-sampX$Ly;Lt<-sampX$Lt
  for (i in 1:n) {
    Ly[[i]]<-(Y[i]-mean(Y))*sampX$Ly[[i]]
  }
  opt<-SetOptions(Ly,Lt,list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,methodSelectK=m1$CV_pi,userSigma2=0.01,usergrid=FALSE))
  obsGrid = sort(unique(c(unlist(Lt))))
  g<-GetSmoothedMeanCurve(Ly,Lt,obsGrid,seq(min(obsGrid), max(obsGrid), length.out = opt$nRegGrid),opt)$muDense
  G<-g%*%Phi_est/l
  b<-G*res$lambda^-1
  beta_pi<-Phi_est%*%t(b)
  res<-FPCA(sampX$Ly,sampX$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m1$CV_In,userSigma2=0.01,usergrid=FALSE))
  Phi_est<-res$phi
  f_In<-lm(Y~res$xiEst)
  beta_In<-Phi_est%*%as.vector(f_In$coef[-1])
  ####1
  res<-FPCA(sampX$Ly,sampX$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m1$CV,userSigma2=0.01,usergrid=FALSE))
  Phi_est<-res$phi
  Xi_1<-matrix(0,nrow = n,ncol = m1$CV)
  for (i in 1:n) {
    for (k in 1:m1$CV) {
      Xi_1[i,k]<-sum(X_grid[i,]*approx(pts_grid,Phi_est[,k],xout = pts_obs[i,],rule = 2)$y)/N
    }
  }
  f1<-lm(Y~Xi_1)
  bVec <- as.vector(f1$coef[-1])
  beta_1<-Phi_est%*%bVec
  out1[1]<-sqrt(trapzRcpp(pts_grid,(beta_1 - beta_true)^2))
  ####pace
  res<-FPCA_CE(sampX$Ly,sampX$Lt,optns =list(nRegGrid=l,methodXi='CE',maxK=M,FVEthreshold=1,userMu=Mu,methodSelectK=m1$CV_pace,userSigma2=0.01,usergrid=FALSE))
  Phi_est<-res$phi
  f_pace<-lm(Y~res$xiEst)
  beta_pace<-Phi_est%*%as.vector(f_pace$coef[-1]) 
  #2
  beta_s1<-matrix(0,nrow = 1,ncol = length(beta_1))
  a1<-c()
  for (s in 1:1) {
    Index_1<-sort(sample.int(n,(n/2)))
    Index_2<-sort(setdiff(c(1:n),Index_1))
    sampX_1<-list(sampX$Lt[Index_1],sampX$Ly[Index_1])
    names(sampX_1)<-c('Lt','Ly')
    sampX_2<-list(sampX$Lt[Index_2],sampX$Ly[Index_2])
    names(sampX_2)<-c('Lt','Ly')
    Y_1<-Y[Index_1]
    Y_2<-Y[Index_2]
    res1<-FPCA(sampX_1$Ly,sampX_1$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m1$CV_split,usergrid=FALSE))
    res2<-FPCA(sampX_2$Ly,sampX_2$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m1$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m1$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m1$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m1$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m1$CV_split) {
        Xi_21[i,k]<-sum(X_grid[Index_1[i],]*approx(pts_grid,Phi_est22[,k],xout = pts_obs[Index_1[i],],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m1$CV_split) {
        Xi_22[i,k]<-sum(X_grid[Index_2[i],]*approx(pts_grid,Phi_est21[,k],xout = pts_obs[Index_2[i],],rule = 2)$y)/N
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
    sampX_1<-list(sampX$Lt[Index_1],sampX$Ly[Index_1])
    names(sampX_1)<-c('Lt','Ly')
    sampX_2<-list(sampX$Lt[Index_2],sampX$Ly[Index_2])
    names(sampX_2)<-c('Lt','Ly')
    Y_1<-Y[Index_1]
    Y_2<-Y[Index_2]
    res1<-FPCA(sampX_1$Ly,sampX_1$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m3$CV_split,usergrid=FALSE))
    res2<-FPCA(sampX_2$Ly,sampX_2$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m3$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m3$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m3$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m3$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m3$CV_split) {
        Xi_21[i,k]<-sum(X_grid[Index_1[i],]*approx(pts_grid,Phi_est22[,k],xout = pts_obs[Index_1[i],],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m3$CV_split) {
        Xi_22[i,k]<-sum(X_grid[Index_2[i],]*approx(pts_grid,Phi_est21[,k],xout = pts_obs[Index_2[i],],rule = 2)$y)/N
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
    sampX_1<-list(sampX$Lt[Index_1],sampX$Ly[Index_1])
    names(sampX_1)<-c('Lt','Ly')
    sampX_2<-list(sampX$Lt[Index_2],sampX$Ly[Index_2])
    names(sampX_2)<-c('Lt','Ly')
    Y_1<-Y[Index_1]
    Y_2<-Y[Index_2]
    res1<-FPCA(sampX_1$Ly,sampX_1$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m5$CV_split,usergrid=FALSE))
    res2<-FPCA(sampX_2$Ly,sampX_2$Lt,optns =list(nRegGrid=l,methodXi='IN',maxK=M,FVEthreshold=1,userMu=Mu,userSigma2=0.01,methodSelectK=m5$CV_split,usergrid=FALSE))
    Phi_est21<-res1$phi
    Phi_est22<-res2$phi
    for (k in 1:m5$CV_split) {
      if(sum((Phi_est21[,k]-Phi_est22[,k])^2)>sum((Phi_est21[,k]+Phi_est22[,k])^2)){Phi_est21[,k]<--Phi_est21[,k]}
    }
    Xi_21<-matrix(0,nrow = n/2,ncol = m5$CV_split)
    Xi_22<-matrix(0,nrow = n/2,ncol = m5$CV_split)
    for (i in 1:(n/2)) {
      for (k in 1:m5$CV_split) {
        Xi_21[i,k]<-sum(X_grid[Index_1[i],]*approx(pts_grid,Phi_est22[,k],xout = pts_obs[Index_1[i],],rule = 2)$y)/N
      }
    }
    for (i in 1:(n/2)) {
      for (k in 1:m5$CV_split) {
        Xi_22[i,k]<-sum(X_grid[Index_2[i],]*approx(pts_grid,Phi_est21[,k],xout = pts_obs[Index_2[i],],rule = 2)$y)/N
      }
    }
    Xi_3<-rbind2(Xi_21,Xi_22)
    Y_3<-c(Y_1,Y_2)
    f3<-lm(Y_3~Xi_3)
    bVec3 <- as.vector(f3$coef[-1])
    beta_s5[s,]<-0.5*(Phi_est21+Phi_est22)%*%bVec3
    a5[s]<-f3$coefficients[1]
  }
  ####
  beta_2<-beta_s1[1,]
  beta_3<-apply(beta_s3, 2, mean)
  beta_5<-apply(beta_s5, 2, mean)
  out1[2]<-sqrt(trapzRcpp(pts_grid,(beta_2 - beta_true)^2))
  out1[3]<-sqrt(trapzRcpp(pts_grid,(beta_3 - beta_true)^2))
  out1[4]<-sqrt(trapzRcpp(pts_grid,(beta_5 - beta_true)^2))
  out1[5]<-sqrt(trapzRcpp(pts_grid,(beta_pi - beta_true)^2))
  out1[6]<-sqrt(trapzRcpp(pts_grid,(beta_pace - beta_true)^2))
  out1[7]<-sqrt(trapzRcpp(pts_grid,(beta_In - beta_true)^2))
  Y_pre1<-c();Y_pre2<-c()
  Y_pre1<-rep(0,2*n);Y_pre2<-rep(0,2*n);Y_pre3<-rep(0,2*n);Y_pre5<-rep(0,2*n);Y_pi<-rep(0,2*n);Y_pace<-rep(0,2*n);Y_In<-rep(0,2*n)
  for (i in 1:(2*n)) {
    temp1<-approx(pts_grid,beta_1,testX$Lt[[i]],rule = 2)$y
    temp2<-approx(pts_grid,beta_2,testX$Lt[[i]],rule = 2)$y
    temp3<-approx(pts_grid,beta_3,testX$Lt[[i]],rule = 2)$y
    tempace<-approx(pts_grid,beta_pace,testX$Lt[[i]],rule = 2)$y
    temp5<-approx(pts_grid,beta_5,testX$Lt[[i]],rule = 2)$y
    tempp<-approx(pts_grid,beta_pi,testX$Lt[[i]],rule = 2)$y
    tempIn<-approx(pts_grid,beta_In,testX$Lt[[i]],rule = 2)$y
    Y_pre1[i]<-f1$coefficients[1]+sum(temp1*testX$Ly[[i]])/N 
    Y_pre2[i]<-a1+sum(temp2*testX$Ly[[i]])/N
    Y_pre3[i]<-mean(a3)+sum(temp3*testX$Ly[[i]])/N
    Y_pre5[i]<-mean(a5)+sum(temp5*testX$Ly[[i]])/N
    Y_pace[i]<-f_pace$coefficients[1]+sum(tempace*testX$Ly[[i]])/N
    Y_In[i]<-f_In$coefficients[1]+sum(tempIn*testX$Ly[[i]])/N
    Y_pi[i]<-sum(tempp*testX$Ly[[i]])/N
  }
  out2[1]<-sum((Y_test-Y_pre1)^2)/(2*n)
  out2[2]<-sum((Y_test -Y_pre2)^2)/(2*n)
  out2[3]<-sum((Y_test -Y_pre3)^2)/(2*n)
  out2[4]<-sum((Y_test -Y_pre5)^2)/(2*n)
  out2[5]<-sum((Y_test -Y_pi)^2)/(2*n)
  out2[6]<-sum((Y_test -Y_pace)^2)/(2*n)
  out2[7]<-sum((Y_test -Y_In)^2)/(2*n)
  out<-list(MSE=out1,Pre=out2,M=m1)
  return(out)
}
#####
system.time({
  cl <- makeCluster(getOption("cl.cores", 50));
  res<-parLapply(cl, 1:200,GetResult)
  stopCluster(cl)
  rm(cl)
})
MSE_10_n100<-matrix(0,nrow = 200,ncol = 7)
Pre_10_n100<-matrix(0,nrow = 200,ncol = 7)
for (r in 1:200) {
  MSE_10_n100[r,]<-res[[r]]$MSE
  Pre_10_n100[r,]<-res[[r]]$Pre
}
apply(MSE_10_n100,2,mean)
apply(Pre_10_n100,2,mean)
apply(MSE_10_n100,2,sd)
apply(Pre_10_n100,2,sd)