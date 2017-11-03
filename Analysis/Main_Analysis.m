
init_agittarget

subjects = [2 3 4 6];

%%
mainTable = zeros(4, 5, length(subjects));
% Columns: 'MeanMS', 'MeanAmp', 'MeanDeltaX', 'MeanDeltaY', 'MeanSaccades'
% Rows: 'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'

for isub = 1:length(subjects)
    tmpTable = table2cell(subjectAnalysis(int2str(subjects(isub))));
    mainTable(:,:,isub) = cell2mat(tmpTable(:,2:6));
end

condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};
meanTable = table(condition_names', zeros(4,1), zeros(4,1), zeros(4,1), 'VariableNames', {'Condition', 'MeanMS', 'MeanAmp', 'MeanSaccades'}); 

for icond = 1:4
    meanTable.MeanMS(icond) = mean(mainTable(icond, 1, :));
    meanTable.MeanAmp(icond) = mean(mainTable(icond, 2, :));
    meanTable.MeanSaccades(icond) = mean(mainTable(icond, 5, :));
end















