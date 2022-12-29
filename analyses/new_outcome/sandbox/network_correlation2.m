function network_correlation2

ss = 1; % which sleep stage (all vs wake vs sleep)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the current output file
data = load([out_folder,'mt_out.mat']);
info = data.out;
all_names = info.all_names;
npts = length(all_names);

all_ic = cell(npts,1);
%% Loop over patients
for p = 1:npts

    name = all_names{p};

    % get the networks
    coh = info.all_coh(p,:,ss);
    re = info.all_re(p,:,ss);
    pc = info.all_pearson(p,:,ss);
    xcor = info.all_xcor(p,:,ss);
    lags = info.all_lags(p,:,ss);
    plv = info.all_plv(p,:,ss);
    nmontages = 3;

    things{1} = coh;
    things{2} = re;
    things{3} = pc;
    things{4} = xcor;
    things{5} = lags;
    things{6} = plv;

    main_net_names = {'coh','re','pc','xcor','lags','plv'};
    montage_names = {'machine','car','bipolar'};

    % Skip if empty
    if isempty(things{1}{1})
        continue; 
    end

    % get number of networks
    nnet = 0;
    for i = 1:length(things)
        nnet = nnet + nmontages*size(things{i}{1},3);
    end

    

    

    % prep inter-network correlation matrix
    ic = nan(nnet,nnet);

    in = 0;
    net_namesi = cell(nnet,1);
    net_namesj = cell(nnet,1);
    for i = 1:length(things)
        for im = 1:nmontages
            for ifreq = 1:size(things{i}{im},3)
                in = in + 1;

                neti = squeeze(things{i}{im}(:,:,ifreq));
                neti = wrap_or_unwrap_adjacency_fc_toolbox(neti);

                net_namesi{in} = sprintf('%s %s %d',main_net_names{i},montage_names{im},ifreq);
                jn = 0;
                % get the jth one
                for j = 1:length(things)
                    for jm = 1:nmontages
                        for jfreq = 1:size(things{j}{jm},3)
                            jn = jn+1;
                            netj = squeeze(things{j}{jm}(:,:,jfreq));
                            netj = wrap_or_unwrap_adjacency_fc_toolbox(netj);
            
                            net_namesj{jn} = sprintf('%s %s %d',main_net_names{j},montage_names{jm},jfreq);
                            ic(in,jn) = corr(neti,netj,'rows','pairwise');
                
                        end
                    end
                end
                
            end
        end
    end

    all_ic{p} = ic;

    if 0
        figure
        turn_nans_gray(ic)
        xticks(1:nnet)
        yticks(1:nnet)
        xticklabels(net_namesj)
        yticklabels(net_namesi)
        colorbar
    end


end

%% Convert to a 3d matrix
% Make all things the biggest size
for i = 1:length(all_ic)
    if isempty(all_ic{i})
        all_ic{i} = nan(size(all_ic{100}));
    end
end

 B=squeeze(cell2mat(cellfun(@(x)reshape(x,[1,1,size(x)]),all_ic,'un',0)));
 avg_mat = squeeze(nanmean(B,1));

figure
turn_nans_gray(avg_mat)
xticks(1:nnet)
yticks(1:nnet)
xticklabels(net_namesj)
yticklabels(net_namesi)
colorbar

end

