
addpath('/net/home/student/l/lgschossmann/git/AgitTarget')

init_agittarget

try
%     newestfile = './cache/data_2017_11_04.mat'; 
    if exist(newestfile,'file')
        tmp = load(newestfile);
        data = tmp.data;
        subjects = str2num(unique(data.subject));
    else
        subjects = [2 3 4 5 6 7 8]; % 
        MS_data = []; % Write in here all microsaccades of all subjects
        SAC_data = [];
        for isub = 1:length(subjects)
            tic
            fprintf('loading %i from %i ...',isub,length(subjects))
            [MStmpTable SACtmpTable] = subjectAnalysis(subjects(isub),1, 1); % the second argument returns a un-aggregated bene-table ;)
            MS_data = [MS_data; MStmpTable];
            SAC_data = [SAC_data; SACtmpTable];
            fprintf(' took %.2f seconds \n',toc)
        end
    end
    
    
%     mainTable = zeros(4, 5, length(subjects));
%     % Columns: 'MeanMS', 'MeanAmp', 'MeanDeltaX', 'MeanDeltaY', 'MeanSaccades'
%     % Rows: 'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'
% 
%     for isub = 1:length(subjects)
%         tmpTable = table2cell(subjectAnalysis(int2str(subjects(isub))));
%         mainTable(:,:,isub) = cell2mat(tmpTable(:,2:6));
%     end
% 
%     condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};
%     meanTable = table(condition_names', zeros(4,1), zeros(4,1), zeros(4,1), 'VariableNames', {'Condition', 'MeanMS', 'MeanAmp', 'MeanSaccades'}); 
% 
%     for icond = 1:4
%         meanTable.MeanMS(icond) = mean(mainTable(icond, 1, :));
%         meanTable.MeanAmp(icond) = mean(mainTable(icond, 2, :));
%         meanTable.MeanSaccades(icond) = mean(mainTable(icond, 5, :));
%     end

catch ME
    if(strcmp(ME.identifier, 'MATLAB:mex:ErrInvalidMEXFile'))
        warning('You need to open matlab via the terminal from the analysis-folder');
    end
end














