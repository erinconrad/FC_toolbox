function show_electrodes_regional_loc


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

%% get variables of interest
locs = data.all_native_locs;
names = data.all_names;
labels = data.all_labels;
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
brainnetome = data.all_brainnetome;
anatomy = data.all_anatomy;
unique_regions = {'L','R'};
nregions = length(unique_regions);


figure
set(gcf,'position',[10 10 1400 800])
tiledlayout(1,3,'tilespacing','tight','padding','tight')
count = 1;
ax_nums = nan(3,1);
cmap = colormap(lines(nregions));
cmap = [0.7 0.7 0.7;cmap];
for p = 71:length(locs)

    

    for ia = 1:3
        switch ia
            case 1 % na atlas approach
                curr_atlas = anatomy{p};
                atlas_lat = label_and_anatomy_lat_determination(labels{p},anatomy{p});
                tit = 'no atlas';
            case 2
                 % AAL approach
                curr_atlas = aal{p};
                atlas_lat = lateralize_regions_simple(aal{p});
                tit = 'aal';
            case 3
                % Brainnetome approach
                curr_atlas = brainnetome{p};
                atlas_lat = lateralize_regions_simple(brainnetome{p});   
                tit = 'brainnetome';
        end
        
        if count == 1
            ax_nums(ia) = nexttile;
        else
            axes(ax_nums(ia))
        end
        hold off
        
        
        [~,ib] = ismember(atlas_lat,unique_regions);

        flocs = locs{p};
        flocs = [flocs;nan nan nan;nan nan nan;nan nan nan];
        ib = [ib;0;1;2];
        curr_atlas = [curr_atlas;{''};{''};{''}];
        curr_labs = labels{p};
        curr_labs = [curr_labs;{''};{''};{''}];
        
        scatter3(flocs(:,1),flocs(:,2),flocs(:,3),150,ib,'filled');
        title(sprintf('Patient %d (%s) atlas %s',p,names{p},tit))
        hold on
        text(flocs(:,1),flocs(:,2),flocs(:,3),...
            cellfun(@(x,y) sprintf('%s (%s)',x,y),curr_atlas,curr_labs,...
            'uniformoutput',false));
        colormap(cmap)
        c = colorbar('ticks',0:nregions,'ticklabels',['none',unique_regions]);
        
    end
    
    pause
    
    count = count+1;
end
    


end