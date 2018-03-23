function errsVals = getAUC_forCV(clssModel, lbls, errsVals, algorithm)
if strcmp(algorithm, 'SVM')
ScoreCVSVMModel = fitSVMPosterior(clssModel);
[~, PostProbs] = kfoldPredict(ScoreCVSVMModel);
else
   [~, PostProbs] =  kfoldPredict(clssModel);
end
[ROCcoord.x,ROCcoord.y,~,AUC_cv] = perfcurve(lbls,PostProbs(:,strcmp(clssModel.ClassNames,'positive')),'positive');
if ROCcoord.x(2) == 0
    ROCcoord.x(2) = [];
    ROCcoord.y(2) = [];
end

errsVals.AUC2 = AUC_cv; errsVals.ROC_coord2 = ROCcoord;
end