QSAR_ridge_regression <- function (X, Y, K) {
     
     # Ridge regression with parcor package
     # Date: 8-4-2014
     # Author: Sheikh M. Rahman
     # require parcor, ggplot2 packages
     # X variables, Y response, k - fold validation
     ridgeFit <- ridge.cv(X, Y, k = K, scale = FALSE); # it will find the optimum lambda by 10-fold validation
     
     coeff_ridge <- matrix(ridgeFit$coefficients,nrow = 1)
     
     
     # fitted data
     yFit <- X %*% t(coeff_ridge) # matrix multiplication
     yFit <- yFit + ridgeFit$intercept # add intercept
     
     residualY <- yFit - Y # residual
     sumResidual <- sum(residualY);
     
     SquaredError <- (yFit - Y)^2
     SumSquareError <- sum(SquaredError);
     #############################################
     # R- square calculation
     # DOF: Model- p (p, no. of explanatory variables, not includes the intercept)
     #  	Error- n-p-1 (n, no. of samples)
     #############################################
     SumSquareModel <- sum((yFit - mean(Y))^2)
     SumSquareTotal <- sum((Y - mean(Y))^2)
     
     Rsquared <- SumSquareModel/SumSquareTotal;
     
     return (list(ridgeFit = ridgeFit, yFit = yFit, residualY = residualY,
                  sumResidual = sumResidual, SquaredError = SquaredError, 
                  SumSquareError = SumSquareError, Rsquared = Rsquared, X = X, Yobs = Y))
     
}


QSAR_lasso_regression <- function (X, Y, K) {
     lassoFit <- adalasso(X,Y, k = K, both = FALSE); # it will find the optimum lambda by k-fold validation
     
     coeff_lasso <- matrix(lassoFit$coefficients.lasso,nrow = 1)
     
     # fitted data
     yFit <- X %*% t(coeff_lasso) # matrix multiplication
     yFit <- yFit + lassoFit$intercept.lasso # add intercept
     yFit[yFit<0] <- 0 # to change the negative value to 0
     residualY <- yFit - Y # residual
     sumResidual <- sum(residualY);
     
     SquaredError <- (yFit - Y)^2
     SumSquareError <- sum(SquaredError);
     
     #############################################
     # R- square calculation
     # DOF: Model- p (p, no. of explanatory variables, not includes the intercept)
     #       Error- n-p-1 (n, no. of samples)
     #############################################
     SumSquareModel <- sum((yFit - mean(Y))^2)
     SumSquareTotal <- sum((Y - mean(Y))^2)
     
     Rsquared <- SumSquareModel/SumSquareTotal;
     
     return (list(lassoFit = lassoFit, yFit = yFit, residualY = residualY,
                  sumResidual = sumResidual, SquaredError = SquaredError, 
                  SumSquareError = SumSquareError, Rsquared = Rsquared, X = X, Yobs = Y))
     
}

plotPanelFigures <- function (Y, regMdl, axsLabel,subHdr) {
yTest_obs <- Y
#coeff <- matrix(regMdl$lassoFit$coefficients.lasso, nrow = 1)
#intercept <- regMdl$lassoFit$intercept.lasso
#yTest <- regMdl$X %*% t( coeff ) # matrix multiplication
yTest <- regMdl$yFit
#yTest[yTest<0] <- 1
#dat <- data.frame(xvar = log(yTest_obs), yvar = log(yTest))
dat <- data.frame(xvar = (yTest_obs), yvar = (yTest))
p <- ggplot(dat, aes(x=xvar, y=yvar)) + labs(title=axsLabel,subtitle=subHdr)+
  geom_point(colour="blue", size = 2.0, shape = 21, fill='orange') +
  theme_bw()  + coord_cartesian(xlim = c(0, 8000), ylim = c(0, 8000)) + #coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  geom_abline(intercept = 0, slope = 1,'linetype' = 'dashed',color = 'black','size' = 0.5) +
  theme (axis.title = element_text(size=21), axis.text=element_text(colour="black",size=10)) +
  labs(x="", y="") +  theme(plot.title = element_text(size=12, hjust=0.5, face="bold", vjust=-1)) +
  theme(plot.subtitle = element_text(size=9, hjust=0.5, vjust=1)) +
  #annotate("text", x=1.25, y=9.5, label= paste("R^2 == ", round(regMdl$Rsquared,2)), parse=TRUE, size = 4,family="Arial", fontface = 2)  
  annotate("text", x=1000, y=7500, label= paste("R^2 == ", round(regMdl$Rsquared,2)), parse=TRUE, size = 4,family="Arial", fontface = "bold")  
  
#  labs(x=paste("Observed", axsLabel), y=paste("Predicted", axsLabel)) 
  
return(p)
}