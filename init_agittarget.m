gitpath = fileparts(which('init_agittarget.m'));
cd(gitpath)

addpath('./lib/edfread/build/linux64');
addpath('./lib/edfread/build/win32');
addpath('./Analysis/')
addpath('./lib/gramm/')

datapath = '/net/store/nbp/projects/agittarget/';
if IsWin && strcmpi(getenv('computername'),'ANANSI')
    % Bene's laptop
    % I use sftp net drive to mount a remote folder on drive "N:"
    datapath = fullfile(['n:' datapath]);

end
