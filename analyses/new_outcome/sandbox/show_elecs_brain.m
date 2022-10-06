function show_elecs_brain

%% parameters
available_hup_ids= [208, 206, 106, 105, 101];

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/new_outcome/data/'];
pt_folder = [locations.main_folder,'data/'];
surf_folder = [locations.main_folder,'data/brain_surf/'];
nishant_code = '~/Desktop/research/general tools/nishant_code/';

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath([nishant_code,'external_lib/']));

%% Load data file
data = load([data_folder,'main_out.mat']);
data = data.out;

%% Load pt file
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;


for sub = 1:length(available_hup_ids)
    id = available_hup_ids(sub);
    
    %% get rid
    for ip = 1:length(pt)
        if strcmp(pt(ip).name,sprintf('HUP%d',id ))
            break
        end
    end
    rid = pt(ip).rid;
    
    %% Get locs
    for jp = 1:length(data.all_names)
        if strcmp(data.all_names{jp},sprintf('HUP%d',id ))
            break
        end
    end
    locs = data.all_native_locs{jp};
    labels = data.all_labels{jp};
   
    subject = ['sub-RID' num2str(rid,'%04.f')];

    [lpv, lpf] = read_surf([surf_folder '/' subject '/surf/lh.pial']);
    [rpv, rpf] = read_surf([surf_folder '/' subject '/surf/rh.pial']);

    %% Figure 1:  Plot pial surface of subject and all electrodes grouped by region
    figure;

    hold on
    hl = trisurf(lpf+1,lpv(:,1),lpv(:,2),lpv(:,3));
    hl.LineStyle = 'none';
    hl.FaceAlpha = 0.1;
    hl.FaceColor = [0.7 0.7 0.7];
    hr = trisurf(rpf+1,rpv(:,1),rpv(:,2),rpv(:,3));
    hr.LineStyle = 'none';
    hr.FaceAlpha = 0.1;
    hr.FaceColor = [0.9 0.7 0.8];
    scatter3(locs(:,1),locs(:,2),...
            locs,100,'w')
    text(locs(:,1),locs(:,2),...
            locs(:,3),labels)

end

end