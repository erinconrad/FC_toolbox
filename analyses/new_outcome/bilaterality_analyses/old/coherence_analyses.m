function hemi_lr = coherence_analyses

%{
Woah bandpower messed up for CAR!!!!! Looks ok for bipolar
I THINK THERE IS AN ERROR IN THIS. SPIKES DONT AGREE, and I think it's pulling out the wrong lateralities 
NEED TO EXCLUDE NON RESECTION
%}

which_atlas = 'aal';
which_outcome = 'ilae';
which_montage = 'car';
which_thing = 'spikes';
symm_cov = 0; % restrict to regions with symmetric coverage
mt_only = 0;
randomize_lats = 0; % set to 1 to randomize soz laterality (null data)

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
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;
soz_lats = data.all_soz_lats;
npts = length(soz_lats);
names = data.all_names;

if strcmp(which_montage,'bipolar') && strcmp(which_thing,'spikes')
    error('why are you doing this')
end
switch which_montage
    case 'bipolar'
        coh = data.all_bipolar_coh;
        fc = data.all_bipolar_fc;
        bp = data.all_bp;
        spikes = data.all_spikes;
    case 'car'  
        coh = data.all_coh;
        fc = data.all_fc;
        bp = data.all_bipolar_bp;
        spikes = data.all_spikes;
end

%% get atlas
switch which_atlas
    case 'aal'
        atlas_names = data.aal_names;
        switch which_montage
            case 'car'
                atlas = data.all_aal;
            case 'bipolar'
                atlas = data.all_bipolar_aal;
        end
    case 'brainnetome'
        atlas_names = data.brainnetome_names;
        switch which_montage
            case 'car'
                atlas = data.all_brainnetome;   
            case 'bipolar'
                atlas = data.all_bipolar_brainnetome;
        end        
end

%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);
outcome = outcome_num;

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Define laterality
old_bilat = strcmp(soz_lats,'bilateral') | strcmp(soz_lats,'diffuse');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');
bilat = nan(length(old_bilat),1);
bilat(old_bilat) = 1;
bilat(unilat) = 0;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

% Define laterality categories
nright = nansum(right_lat); nleft = nansum(left_lat); nbilat = nansum(bilat); nempty = sum(cellfun(@isempty,soz_lats));
assert(nright+nleft+nbilat+nempty == npts)
lat_num = nan(npts,1);
lat_num(cellfun(@isempty,soz_lats)) = 3;
lat_num(bilat==1) = 0;
lat_num(left_lat==1) = 1;
lat_num(right_lat==1) = 2;
assert(~any(isnan(lat_num)))

%% Randomize lats (null data)
% My coherence results do indeed turn non significant for the SOZ vs non
% SOZ laterality side
if randomize_lats
    new_lat_num = lat_num(randperm(npts));
    bilat = nan(length(old_bilat),1);
    bilat(new_lat_num == 0) = 1;
    bilat(ismember(new_lat_num,[1 2])) = 0;

    right_lat = zeros(npts,1);
    right_lat(new_lat_num==2) = 1;

    left_lat = zeros(npts,1);
    left_lat(new_lat_num==1) = 1;

    assert(nansum(bilat)==nbilat && sum(right_lat)==nright && sum(left_lat)==nleft)

end

%% Decide thing
switch which_thing
    case 'fc'
        thing = fc;
    case 'coh'
        thing = coh;
    case 'bp'
        thing = bp;
    case 'spikes'
        thing = spikes;
    case 'nelecs'
end

%% Put into atlas
thing_atlas = convert_to_atlas(thing,atlas,atlas_names);
dims = ndims(thing_atlas);

%% Error checking - do atlas parcellations agree
for ip = 1:npts
    switch dims
        case 2
            non_nan_regions = (~isnan(squeeze(thing_atlas(ip,:))))';
        case 3
            non_nan_regions = (~isnan(squeeze(nanmean(thing_atlas(ip,:,:),3))))';
        case 4
            non_nan_regions = (~isnan(squeeze(nanmean(thing_atlas(ip,:,:,:),[3 4]))))';
    end
    if ~(isequal(unique(atlas_names(non_nan_regions)),unique(atlas{ip}(cellfun(@(x) ~isempty(x),atlas{ip})))))
        fprintf('\nsuprise nans for %s\n',names{ip});
    end
    %setdiff(unique(atlas_names(non_nan_regions)),unique(atlas{ip}(cellfun(@(x) ~isempty(x),atlas{ip}))))
end

%% Atlas lateralizations
[locs,lats] = lateralize_regions(atlas_names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither_lat = ~left & ~right;

% confirm atlas has as many on right as left
assert(sum(left)==sum(right));

%% Re-order atlas to be left then right then neither
lr_order = reorder_lr(locs,lats); % get the order

% Put everything into this order
left = left(lr_order);
right = right(lr_order);
neither_lat = neither_lat(lr_order);
atlas_names = atlas_names(lr_order);
switch dims
    case 2
        thing_atlas = thing_atlas(:,lr_order);
    case 3
        if size(thing_atlas,2) == size(thing_atlas,3)
            thing_atlas = thing_atlas(:,lr_order,lr_order);
        else
            thing_atlas = thing_atlas(:,lr_order,:);
        end
    case 4
        thing_atlas = thing_atlas(:,lr_order,lr_order,:);
end
locs = locs(lr_order);
lats = lats(lr_order);


%% index of contralateral region for each region
contra_index = nan(size(left));
contra_index(1:sum(left)) = ([1:sum(left)])'+sum(left); % the right sided ones start after left, so contralateral ones for the first nleft are these
contra_index(sum(left)+1:sum(left)*2) = ([1:sum(left)])';

% make sure first half locs are same as last half locs
assert(isequal(locs(1:sum(left)),locs(sum(left)+1:sum(left)*2)))
assert(isequal(locs(contra_index(1:sum(left))),locs(contra_index(sum(left)+1:sum(left)*2))))

%% Build symmetric coverage atlas
if symm_cov
    [thing_atlas,all_bilateral] = symmetric_coverage_atlas(thing_atlas,locs,lats);
end

%% Restrict to mt
% Define names corresponding to mesial temporal

if mt_only
    switch which_atlas
        case 'aal'
            mt_names = {'Hippocampus','Amygdala'};
        case 'brainnetome'
            mt_names = {'Amyg','Hipp'};
    end
    
    mt = contains(atlas_names,mt_names);
    
    switch dims
        case 2
            thing_atlas(:,~mt) = nan;
        case 3
            if size(thing_atlas,2) == size(thing_atlas,3)
                thing_atlas(:,~mt,:) = nan; thing_atlas(:,:,~mt) = nan;
            else
                thing_atlas(:,~mt,:) = nan;
            end
        case 4
            thing_atlas(:,~mt,:,:) = nan; thing_atlas(:,:,~mt,:) = nan;
    end

end
    
   

%% Get average left and right
switch dims
    case 2
        hemi_lr = nan(npts,2,1);
    case 3
        if size(thing_atlas,2) == size(thing_atlas,3)
            hemi_lr = nan(npts,2,1);
        else
            hemi_lr = nan(npts,2,size(thing_atlas,3));
        end
    case 4
        hemi_lr = nan(npts,2,size(thing_atlas,4));
end
for ip = 1:npts
    switch dims
        case 2
            hemi_lr(ip,:,:) = [nanmean(thing_atlas(ip,left),'all'),...
                nanmean(thing_atlas(ip,right),'all')];
        case 3
            if size(thing_atlas,2) == size(thing_atlas,3)
                hemi_lr(ip,:,:) = [nanmean(thing_atlas(ip,left,left),'all'),...
                    nanmean(thing_atlas(ip,right,right),'all')];
            else
                hemi_lr(ip,:,:) = [nanmean(thing_atlas(ip,left,:),[1 2]),...
                    nanmean(thing_atlas(ip,right,:),[1 2])];
            end
        case 4
            hemi_lr(ip,:,:) = [nanmean(thing_atlas(ip,left,left,:),[1 2 3]),...
                nanmean(thing_atlas(ip,right,right,:),[1 2 3])];
    end

end

return

%% Get average SOZ side - non SOZ side
hemi_soz_non = nan(size(hemi_lr));
for ip = 1:npts
    if left_lat(ip) % keep same side
        switch dims
            case 2
                hemi_soz_non(ip,:,:) = hemi_lr(ip,:,:);
            case 3
                hemi_soz_non(ip,:,:) = hemi_lr(ip,:,:);
            case 4
                hemi_soz_non(ip,:,:) = hemi_lr(ip,:,:);
        end
    elseif right_lat(ip) % swap 1st and 2nd element in 2nd dimension
        switch dims
            case 2
                hemi_soz_non(ip,:,:) = [hemi_lr(ip,2,:) hemi_lr(ip,1,:)];
            case 3
                hemi_soz_non(ip,:,:) = [hemi_lr(ip,2,:) hemi_lr(ip,1,:)];
            case 4
                hemi_soz_non(ip,:,:) = [hemi_lr(ip,2,:) hemi_lr(ip,1,:)];
        end

    end


end

%% Compare intra-hemispheric thing on SOZ vs non SOZ side
% For Pearson correlation, this gives same result as in JNE paper
if 0
    figure
    for f = 1:size(hemi_soz_non,3)
        nexttile
        paired_plot(hemi_soz_non(:,:,f),'thing',{'soz','non-soz'})
    end
    
end

%% Calculate asymmetry index of L-R
ai = squeeze(asymmetry_index(hemi_lr(:,1,:),hemi_lr(:,2,:)));

%% Compare asymmetry index between those with unilateral vs bilateral SOZ
if 0
    figure
    for f = 1:size(ai,2)
        nexttile
        unpaired_plot(ai(unilat==1,f),ai(bilat==1,f),{'unilateral','bilateral'},'thing');
    end
end

%% Compare asymmetry index between those with good vs bad outcome
if 1
    figure; set(gcf,'position',[100 100 300*size(ai,2) 350])
    tiledlayout(1,size(ai,2))
    for f = 1:size(ai,2)
        nexttile
        unpaired_plot(ai(outcome==1,f),ai(outcome==0,f),{'good outcome','bad outcome'},'thing');
    end
end

end

