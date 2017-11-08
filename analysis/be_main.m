
init_agittarget


newestfile = './cache/data_2017_11_04.mat';
if exist(newestfile,'file')
    tmp = load(newestfile);
    data = tmp.data;
else
    subjects = [2 3 4 6];
    data = [];
    for isub = 1:length(subjects)
        tic
        fprintf('loading %i from %i ...',isub,length(subjects))
        tmpTable = subjectAnalysis(subjects(isub),1); % the second argument returns a un-aggregated bene-table ;)
        data = [data; tmpTable];
        fprintf(' took %.2f seconds \n',toc)
    end
end

%%
data.subject = str2num(data.subject);
data(data.Amplitude>10^3,:) = [];
%%
figure
g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
g.geom_point();
g.facet_grid([],data.condition)
g.draw()
%%
figure
g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
g.stat_smooth();
g.draw()

%% Now let's have a look at more interesting things
% amplitude densities
figure
g = gramm('x',log10(data.Amplitude),'color',data.condition);

g.facet_grid(data.subject,[]);
g.stat_density();
g.draw();
% set(g.facet_axes_handles,'XLim',[0 70])

xlim([0 70]);

%% number of MS per trial
stat = grpstats(data,{'subject','trial','condition'});
figure
g = gramm('x',stat.condition,'y',stat.GroupCount,'color',stat.condition);

g.stat_violin();
g.geom_jitter();
g.facet_grid([],stat.subject);
g.stat_summary('geom','black_errorbar');
g.draw();
set(findobj(g.facet_axes_handles,'Marker','o'),'MarkerFaceColor',[.2 .2 .2])