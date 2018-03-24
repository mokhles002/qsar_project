###################################################################
# run functions of lasso regression to predict continous endpoint
###################################################################
# import the require packages
library("parcor")
library("ggplot2")
library("ggpubr")
# import the necessary functions for regression and plotting
source("lib/QSAR_functions.r")
# import the data
load("data/EColi_TELI_Physical_EC50_23chem_5pathways.rdata")

# lasso regression from function
Y <- EC50               # Response variable
X <-  X_TELI_total      # Explanatory variables
descriptorType <- "total"     # filename suffix to save the results
# X <- cbind(X_TELI_total,X_physical) # To generate explanatory variables for combined model
# descriptorType <- "physical_total"
lasso_regression_summary <- QSAR_lasso_regression (X, Y, 10)
"R-square : "; lasso_regression_summary$Rsquared

# save the results as rdata
save(file = paste("ec50_",descriptorType,"_Rsquare_",lasso_regression_summary$Rsquared,".rdata",sep = ''), lasso_regression_summary,X,Y)

#########################################################
# Draw figure for 1 model
#########################################################

# prepare data for plotting from raw data and prediction models
xTest <- X
yTest_obs <- Y
coeff <- matrix(lasso_regression_summary$lassoFit$coefficients.lasso, nrow = 1)
intercept <- lasso_regression_summary$lassoFit$intercept.lasso

yTest <- xTest %*% t( coeff ) # matrix multiplication
yTest <- yTest + intercept # add intercept

dat <- data.frame(xvar = log(yTest_obs), yvar = log(yTest))

# save the observed vs. predicted for one type of explanatory variables
tiff(file = paste("vF_EC50_",descriptorType,".tiff",sep = ''), width = 5.0, height = 4.5, units = 'in', res = 300, compression = 'lzw')
ggplot(dat, aes(x=xvar, y=yvar)) + geom_point(colour="blue", size = 4, shape=1) +
  theme_bw()  + coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  geom_abline(intercept = 0, slope = 1,'linetype' = 'dashed',color = 'black','size' = 0.5) +
  theme (axis.title = element_text(size=21), axis.text=element_text(colour="black",size=17)) + 
  labs(x="Observed EC50", y="Predicted EC50") 
dev.off()


############################################################################
###     Plot figures of all 7 models in a panel
############################################################################

# load the results of the 7 models
load ("data/ec50_physical_Rsquare_0.130349046857749.rdata")
physicalModel = lasso_regression_summary
physicalModel$X = X
load ("data/ec50_gene_Rsquare_0.704734052060167.rdata")
geneModel = lasso_regression_summary
geneModel$X = X
load ("data/ec50_pathways_Rsquare_0.783956863877711.rdata")
pathModel = lasso_regression_summary
pathModel$X = X
load ("data/ec50_total_Rsquare_0.0848274215748496.rdata")
totalModel = lasso_regression_summary
totalModel$X = X
load("data/ec50_physical_gene_Rsquare_0.318469473082123.rdata")
genePhysicalModel = lasso_regression_summary
genePhysicalModel$X = X
load("data/ec50_physical_pathways_Rsquare_0.130349046857749.rdata")
pathPhysicalModel = lasso_regression_summary
pathPhysicalModel$X = X
load("data/ec50_physical_total_Rsquare_0.130349046857749.rdata")
totalPhysicalModel = lasso_regression_summary
totalPhysicalModel$X = X

# generate the panel plots for all 7 models

p1 <- plotPanelFigures(Y,physicalModel,'QSAR model','(Physio-chemical descriptors)') + labs(x=NULL, y =NULL)
p2 <- plotPanelFigures(Y,geneModel,'QBAR-Biomarker model','(Biomarker level descriptors)') + labs(x=NULL, y =NULL)
p3 <- plotPanelFigures(Y,genePhysicalModel,'QSBAR-Biomarker model','(Physio-chemical & biomarker level descriptors)') + labs(x=NULL, y =NULL)
p4 <- plotPanelFigures(Y,pathModel,'QBAR-Pathway model','(Pathway level descriptors)') + labs(x=NULL, y =NULL)
p5 <- plotPanelFigures(Y,pathPhysicalModel,'QSBAR-Pathway model','(Physio-chemical & pathway level descriptors)') + labs(x=NULL, y =NULL)
p6 <- plotPanelFigures(Y,totalModel,'QBAR-Cellular model','(Cellular level descriptors)') + labs(x=NULL, y =NULL)
p7 <- plotPanelFigures(Y,totalPhysicalModel,'QSBAR-Cellular model','(Physio-chemical & cellular level descriptors)') + labs(x=NULL, y =NULL)

# save the panel figures in the same image file as "tiff"
tiff(file = "allFigures_combined_5pathways_nolog_vF.tiff", width = 6.5, height = 8, units = 'in', res = 300, compression = 'lzw')
fig <- ggarrange(p1, p2,p4,p6, p3, p5, p7,
          labels = c("a)", "b)", "c)","d)", "e)", "f)", "g)"),
          ncol = 2, nrow = 4,hjust = -1)

annotate_figure(fig,
                bottom = text_grob(expression(bold("Observed EC"[50])), color = "black",
                                    face = "bold", size = 14),
                left = text_grob(expression(bold("Predicted EC"[50])), color = "black", 
                                 face = "bold",rot = 90, size = 14)
                )
dev.off()
##############################################################################
#   Save the model results and coefficients as csv file.
##############################################################################

regMdl <- lasso_regression_summary
coeff_lasso <- matrix(regMdl$lassoFit$coefficients.lasso,nrow = 1)
path_coeff <- matrix(coeff_lasso[,1:42], nrow = 6, byrow = TRUE)# reorder the coefficient as conc. X pathways
write.csv(path_coeff, file = "lasso_coefficient_path_physical_onlyPath.csv")
write.csv(coeff_lasso[,-(1:42)], file = "lasso_coefficient_path_physical_onlyPhysical.csv")
write.csv(t(coeff_lasso), file = "lasso_coefficient_physical.csv")
write.csv(physicalDescriptorName, file = "physical_descriptorName.csv")
write.csv(t(geneName), file = "gene_descriptorName.csv")


