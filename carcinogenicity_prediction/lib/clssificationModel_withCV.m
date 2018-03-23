function [errsVals, clssModel] = clssificationModel_withCV(datavalues, lbls, algorithm,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform the gene selection from the genotoxicity
%   Developed by: Sheikh M. Rahman
%   Date: 02-25-2016
%   A fucntion that find the error using SVM classifier
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_pos = strcmp(lbls,'positive');
id_neg = strcmp(lbls,'negative');

dataTrainNegative = datavalues(id_neg,:);
dataTrainPositive = datavalues(id_pos,:);
lblNegative = lbls(id_neg,:);
lblPositive = lbls(id_pos,:);

switch algorithm
    case 'knn_3'
        clssModel = fitcknn (datavalues, lbls, 'NumNeighbors',3,'kfold', k);
    case 'knn_5'
        clssModel = fitcknn (datavalues, lbls,'kfold', k, 'NumNeighbors',5);
    case 'SVM'
        clssModel = fitcsvm(datavalues,lbls,'kfold', k);
end



yP = kfoldPredict(clssModel);
accuracy = 1-kfoldLoss(clssModel);
sensitivity = sum(strcmp(lblPositive,yP(id_pos)))/sum(id_pos);
specificity = sum(strcmp(lblNegative,yP(id_neg)))/sum(id_neg);
errsVals.accuracy = accuracy; errsVals.sensitivity = sensitivity; errsVals.specificity = specificity;

% AUC

if strcmp(algorithm, 'SVM')
ScoreCVSVMModel = fitSVMPosterior(clssModel);
[~, PostProbs] = kfoldPredict(ScoreCVSVMModel);
else
   [~, PostProbs] =  kfoldPredict(clssModel);
end
[ROCcoord.x,ROCcoord.y,~,AUC_cv] = perfcurve(lbls,PostProbs(:,strcmp(clssModel.ClassNames,'positive')),'positive', 'XVals',0:0.05:1);
if ROCcoord.x(2) == 0
    ROCcoord.x(2) = [];
    ROCcoord.y(2) = [];
end

errsVals.AUC = AUC_cv; errsVals.ROC_coord = ROCcoord;
end
