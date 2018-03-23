function plotAUC_function (errsVals, clr)
%plot(errsVals.ROC_coord.x,errsVals.ROC_coord.y,'color',clr,'linestyle','--')
hold on;
plot(errsVals.ROC_coord2.x,errsVals.ROC_coord2.y,'color',clr,'linestyle','-')
scatter(1-errsVals.specificity,errsVals.sensitivity,'markeredgecolor',clr);
end