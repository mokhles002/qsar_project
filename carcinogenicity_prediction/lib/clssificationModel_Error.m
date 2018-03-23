function [errsVals, clssModel] = clssificationModel_Error(datavalues, lbls, algorithm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform the gene selection from the genotoxicity
%   Developed by: Sheikh M. Rahman
%   Date: 02-25-2016
%   A fucntion that find the error using SVM classifier
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_positive = strcmp(lbls,'positive');
id_negative = strcmp(lbls,'negative');

dataTrainNegative = datavalues(id_negative,:);
dataTrainPositive = datavalues(id_positive,:);
lblNegative = lbls(id_negative,:);
lblPositive = lbls(id_positive,:);

switch algorithm
    case 'knn_3'
        clssModel = fitcknn (datavalues, lbls, 'NumNeighbors',3);
    case 'knn_5'
        clssModel = fitcknn (datavalues, lbls, 'NumNeighbors',5);
    case 'SVM'
        clssModel = fitcsvm(datavalues,lbls);
end

accuracy = (1-loss(clssModel, datavalues,lbls))*100;
sensitivity = (1-loss(clssModel, dataTrainPositive, lblPositive))*100;
specificity =(1-loss(clssModel, dataTrainNegative, lblNegative))*100;

errsVals.accuracy = accuracy; errsVals.sensitivity = sensitivity; errsVals.specificity = specificity;
end
