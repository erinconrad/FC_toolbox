function info_on_designations

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];


% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

out_folder1 = [scripts_folder,'analyses/sleep/data/'];

% Load out folder
out = load([out_folder1,'out.mat']);
out = out.out;

swdes = out.roc_out.swdes;
nawake = 0;
nasleep = 0;
nuncertain = 0;
nn1 = 0;
nn2 = 0;
nr = 0;
nn3 = 0;
na = 0;

for i = 1:length(swdes)
    nawake = nawake + length(swdes(i).sw.wake);
    nasleep = nasleep + length(swdes(i).sw.sleep);
    nuncertain = nuncertain + length(swdes(i).sw.uncertain);

    na = na + sum(strcmp(swdes(i).designations,'a'));
    nn1 = nn1 + sum(strcmp(swdes(i).designations,'1'));
    nn2 = nn2 + sum(strcmp(swdes(i).designations,'2'));
    nn3 = nn3 + sum(strcmp(swdes(i).designations,'3'));
    nr = nr + sum(strcmp(swdes(i).designations,'r'));
end

end