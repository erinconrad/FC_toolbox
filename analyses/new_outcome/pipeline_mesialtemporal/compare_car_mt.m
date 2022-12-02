function [all_pc_corr,all_coh_corr] = compare_car_mt

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
edf_path = [results_folder,'edf_out/'];

%% allowable labels
allowable_labels = get_allowable_elecs;

%% Load the main data file
main = load([inter_folder,'main_out.mat']);
main = main.out;
names = main.all_names;
npts = length(names);

all_pc_corr = nan(npts,1);
all_coh_corr = nan(npts,6);

% Loop over pts
for ip = 1:npts
    name = names{ip};

    % Load the corresponding edf summary
    if exist([edf_path,name,'/summ.mat'],'file') ==0
        continue
    end
    edf_summ = load([edf_path,name,'/summ.mat']);
    edf_summ = edf_summ.out;

    % Get the allowed labels for the main thing
    labels = main.all_labels{ip,1};
    allowed = ismember(labels,allowable_labels);
    main_labels = labels(allowed);

    % get the labels from the edf thing
    edf_labels = edf_summ.labels;

    % make sure theyre the same (they might not be if they added MT elecs
    % later)
    if ~(isequal(main_labels,edf_labels))
        fprintf('\nUnequal elecs for %s\n',name);
        main_labels
        edf_labels
        continue
    end

    % reduce main to allowed
    coh = main.all_coh{ip,1,1};
    pc = main.all_pearson{ip,1,1};
    coh = coh(allowed,allowed,:);
    pc = pc(allowed,allowed);

    edf_pc = edf_summ.all_pc;
    edf_coh = edf_summ.all_coh;
    edf_pc = squeeze(nanmean(edf_pc,1));
    edf_coh = squeeze(nanmean(edf_coh,1));

    % wrap for the purpose of correlating
    pc_wrap = wrap_or_unwrap_adjacency_fc_toolbox(pc);
    coh_wrap = wrap_or_unwrap_adjacency_fc_toolbox(coh);
    edf_pc_wrap = wrap_or_unwrap_adjacency_fc_toolbox(edf_pc);
    edf_coh_wrap = wrap_or_unwrap_adjacency_fc_toolbox(edf_coh);

    % correlate
    all_pc_corr(ip) = corr(pc_wrap,edf_pc_wrap ,'rows','pairwise');
    A = corr(coh_wrap,edf_coh_wrap ,'rows','pairwise');
    all_coh_corr(ip,:) = diag(A);

    if 0
        figure
        nexttile
        turn_nans_gray(pc)
        nexttile
        turn_nans_gray(edf_pc)
    end
end

end