function compare_atlases

%% Parameters
atlases = {'aal','aal_bernabei','aal_bernabei_bipolar'};


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Prep
natlases = length(atlases);
atlas_names_ns = cell(natlases,2);

for in = 1:natlases
    out = show_atlas(atlases{in},0);
    ns = out.ns;
    names = out.names;
    atlas_names_nums{in,2} = ns;
    atlas_names_nums{in,1} = names;
    
end

% Build pairs
npairs = nchoosek(natlases,2);
pairs = nan(npairs,2);
which_pair = 0;
for in = 1:natlases
    for im = in+1:natlases
        which_pair = which_pair+1;
        pairs(which_pair,:) = [in im];
    end
end

% Correlate
for ip = 1:npairs
    
   
    labels1 = atlas_names_nums{pairs(ip,1),1};
    labels2 = atlas_names_nums{pairs(ip,2),1};
    ns1 = atlas_names_nums{pairs(ip,1),2};
    ns2 = atlas_names_nums{pairs(ip,2),2};

    [labels1,ns1,labels2,ns2] = find_matching_labels(labels1,ns1,labels2,ns2);
    r = corr(ns1,ns2);
     
    fprintf(['\nComparing %s and %s:'...
        ' ns r = %1.2f\n'],atlases{pairs(ip,1)},atlases{pairs(ip,2)},...
        r);
    
    
end

end