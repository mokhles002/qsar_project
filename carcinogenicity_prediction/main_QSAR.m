%% Prediction of the in-vivo carcinogenicity endpoint
% the outcome is predicted using the SVM classification algorithm
% data_descriptor.mat contains the physio-chemical properties as well as
% the biological descriptors

load ('data/data_descriptor.mat'); % load the pre-stored data.
addpath lib;
clssAlgorithm = 'SVM';
kfld = 10; %n=10;
%rng(n,'twister');
[errsVals_physical, clssModel_physical] = clssificationModel_withCV(physical_descriptor, clsslbl, clssAlgorithm, kfld);
[errsVals_gene, clssModel_gene] = clssificationModel_withCV(gene_descriptor, clsslbl, clssAlgorithm, kfld);
[errsVals_path, clssModel_path] = clssificationModel_withCV(path_descriptor, clsslbl, clssAlgorithm, kfld);
[errsVals_total, clssModel_total] = clssificationModel_withCV(total_descriptor, clsslbl, clssAlgorithm, kfld);
[errsVals_gene_physical, clssModel_gene_physical] = clssificationModel_withCV([gene_descriptor, physical_descriptor], clsslbl, clssAlgorithm, kfld);
[errsVals_path_physical, clssModel_path_physical] = clssificationModel_withCV([path_descriptor, physical_descriptor], clsslbl, clssAlgorithm, kfld);
[errsVals_total_physical, clssModel_total_physical] = clssificationModel_withCV([total_descriptor, physical_descriptor], clsslbl, clssAlgorithm, kfld);

%{
errsVals_physical
errsVals_gene
errsVals_path
errsVals_total
errsVals_gene_physical
errsVals_path_physical
errsVals_total_physical
%}
%
save('results/svm_carcinogen_vF_CV.mat','errsVals_gene','errsVals_path',...
    'errsVals_physical','errsVals_total','clssModel_gene', 'clssModel_path',...
    'clssModel_physical', 'clssModel_total','errsVals_gene_physical',...
    'errsVals_path_physical','errsVals_total_physical','clssModel_gene_physical',...
    'clssModel_path_physical','clssModel_total_physical');
%% Plot Heatmap of the prediction label

% get data for the heatmap
% find the outcome from each prediction model 
predictedLabel = [kfoldPredict(clssModel_physical),...
    kfoldPredict(clssModel_gene),...
    kfoldPredict(clssModel_path),...
    kfoldPredict(clssModel_total),...
    kfoldPredict(clssModel_gene_physical),...
    kfoldPredict(clssModel_path_physical),...
    kfoldPredict(clssModel_total_physical),...
    ];

% determine whether the outcome label is true or false
for i = 1:7
    labelTrueFalse(:,i) = strcmp(predictedLabel(:,i),clsslbl);
end
mthd = {'QSAR'; 'QBAR-Biomarker'; 'QBAR-Pathway'; 'QBAR-Cellular';...
            'QSBAR-Biomarker'; 'QSBAR-Pathway'; 'QSBAR-Cellular';};
save ('labelForHeatmap.mat','predictedLabel','clsslbl','chemName', 'labelTrueFalse','mthd');

% generate the  Heatmap

id_positive = strcmp(clsslbl,'positive');
id_negative = strcmp(clsslbl,'negative');
nPos = sum(id_positive); nNeg = sum(id_negative);
labelTrueFalse_ord = [labelTrueFalse(id_positive,:); labelTrueFalse(id_negative,:)];
chemName_ord = [chemName(id_positive); chemName(id_negative)];
fh = figure;
set(fh, 'PaperUnits','inches', 'PaperSize',[6.5 7], 'PaperPosition',[0 0 6.5 7],...
    'Units','inches', 'Position',[0 0 6.5 7],'Color', 'w');
hmain = axes('position', [0.36 0.315 0.62 0.68]);
imagesc(labelTrueFalse_ord);
colormap([[.85 0 0];[0 0.85 0]]);

set(gca,'xtick',0.5:1:7.5, 'xticklabel','','xticklabelrotation',90,...
    'ytick',0.5:1:20.5, 'yticklabel','', 'tickdir','out', 'fontname','arial','fontsize', 12);
set(gca, 'ygrid','on','xgrid','on','layer', 'top','GridColor','w','gridalpha',1,...
    'gridlinestyle','--', 'linewidth',0.7);
set(gca, 'ticklength', [0.005 0.1]);

t1 = text(repmat(0.38,nPos,1), 1:1:nPos, chemName_ord(1:nPos),'color','b',...
    'fontname', 'arial', 'fontsize', 12);
set(t1, 'HorizontalAlignment','right','VerticalAlignment','middle')
t2 = text(repmat(0.38,nNeg,1), nPos+1:1:nPos+nNeg, chemName_ord(nPos+1:1:end),...
    'fontname', 'arial', 'fontsize', 12);
set(t2, 'HorizontalAlignment','right','VerticalAlignment','middle')

t2 = text(1:1:7,repmat(20.7,7,1), mthd,...
    'fontname', 'arial', 'fontsize', 13);
set(t2, 'HorizontalAlignment','right','VerticalAlignment','middle','rotation',90)

annotation('textbox',[0.37 0.05 0.62 0.03],'str','Prediction Model Type',...
            'linestyle','none','margin',0,'fontname','Arial','fontWeight','bold',...
            'fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle');
ha = axes('position',[0.018 0.502 0.05 0.492]);
drawbrace([0 2],[0 20],0.2,'color','b');ylim([0 20]);ha.Visible = 'off';
text(-2.15, 9.2, 'Positive','rotation',90,...
    'fontname','Arial','fontsize',15, 'fontWeight','bold');

hb = axes('position',[0.018 0.295 0.05 0.25]);
drawbrace([0 2],[0 20],0.2,'color','k');ylim([0 20]);hb.Visible = 'off';
text(-2.15, 7.5, 'Negative','rotation',90,...
    'fontname','Arial','fontsize',15, 'fontWeight','bold');

ha = axes('position',[0.37 0.03 0.62 0.1]);
hb=bar([1:5;1:5]', 'edgecolor','none'); colormap(ha,[[0 0.85 0];[.85 0 0]]);
gridLegend(hb,2,{'Correct prediction', 'Wrong prediction'},'box','off',...
    'fontname','Arial','fontsize',12,'location','northoutside');
cla(ha); set(gca, 'Visible', 'off');
set(gcf, 'InvertHardCopy', 'off'); % setting 'grid color reset' off
print(fh,'-dpdf','-r300','results/heatmap_carcinogen_prediction_vF.pdf')

%% Export the prediction performance parameters as CSV format
accrcy_all = [errsVals_physical.accuracy; errsVals_gene.accuracy;
    errsVals_path.accuracy; errsVals_total.accuracy; errsVals_gene_physical.accuracy;
    errsVals_path_physical.accuracy; errsVals_total_physical.accuracy]*100;

snstivty_all = [errsVals_physical.sensitivity; errsVals_gene.sensitivity;
    errsVals_path.sensitivity; errsVals_total.sensitivity; errsVals_gene_physical.sensitivity;
    errsVals_path_physical.sensitivity; errsVals_total_physical.sensitivity;]*100;

spfcty_all = [errsVals_physical.specificity; errsVals_gene.specificity;
    errsVals_path.specificity; errsVals_total.specificity; errsVals_gene_physical.specificity;
    errsVals_path_physical.specificity; errsVals_total_physical.specificity;]*100;

AUC_all = [errsVals_physical.AUC; errsVals_gene.AUC;
    errsVals_path.AUC; errsVals_total.AUC; errsVals_gene_physical.AUC;
    errsVals_path_physical.AUC; errsVals_total_physical.AUC;]*100;

tableForExport = [array2table(mthd,...
'VariableNames',{'model_type'}), array2table([AUC_all/100,  accrcy_all, snstivty_all, spfcty_all],...
    'variableNames',{'AUC','Accuracy','Sensitivity','Specificty'})];
writetable(tableForExport,'results/summary_with_svm.csv');





