function network_correlation

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_summ_out/'];
sleep_stage_path = [results_folder,'edf_out/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the current output file
data = load([out_folder,'main_out.mat']);
data = data.out;
all_names = data.all_names;
npts = length(all_names);

all_ic = cell(npts,1);
%% Loop over patients
for p = 1:npts

    name = all_names{p};

    % Load the edf summary file
    if exist([edf_path,name,'/summ.mat'],'file') == 0, continue; end
    info = load([edf_path,name,'/summ.mat']);
    info = info.out;

    % get the networks
    coh = info.all_coh;
    re = info.all_re;
    pc = info.all_pc;
    xcor = info.all_xcor;
    lags = info.all_lags;
    plv = info.all_plv;
    nmontages = 3;

    things{1} = coh;
    things{2} = re;
    things{3} = pc;
    things{4} = xcor;
    things{5} = lags;
    things{6} = plv;

    main_net_names = {'coh','re','pc','xcor','lags','plv'};
    
    % average across time
    for i = 1:length(things)
        things{i} = squeeze(nanmean(things{i},1));
    end

    % get number of networks
    nnet = 0;
    for i = 1:length(things)
        nnet = nnet + size(things{i},1)*size(things{i},4);
    end

    % prep inter-network correlation matrix
    ic = nan(nnet,nnet);
    in = 0;
    net_namesi = cell(nnet,1);
    net_namesj = cell(nnet,1);
    for i = 1:length(things)
        for im = 1:nmontages
            for ifreq = 1:size(things{i},4)
                in = in + 1;

                neti = squeeze(things{i}(im,:,:,ifreq));
                neti = wrap_or_unwrap_adjacency_fc_toolbox(neti);

                net_namesi{in} = sprintf('%s %s %d',main_net_names{i},info.montages{im,1},ifreq);
                jn = 0;
                % get the jth one
                for j = 1:length(things)
                    for jm = 1:nmontages
                        for jfreq = 1:size(things{j},4)
                            jn = jn+1;
                            netj = squeeze(things{j}(jm,:,:,jfreq));
                            netj = wrap_or_unwrap_adjacency_fc_toolbox(netj);
            
                            net_namesj{jn} = sprintf('%s %s %d',main_net_names{j},info.montages{jm,1},jfreq);
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

all_ic = cell2mat(all_ic);

end

