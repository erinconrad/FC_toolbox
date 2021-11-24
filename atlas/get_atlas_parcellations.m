function out = get_atlas_parcellations(rid,elabels,name)

%% Parameters
which_atlas = 'aal';

%% Get file locs
locations = fc_toolbox_locs;
atlas_folder = [locations.main_folder,'data/atlas/'];
parcellation_folder = [atlas_folder,'atlas_localizations/'];

%% Clean labels
elabels = decompose_labels(elabels,name);

%% Load correct atlas table
if rid<100
    ridtext = sprintf('00%d',rid);
else
    ridtext = sprintf('0%d',rid);
end

listing = dir(sprintf('%s*%s*.csv',parcellation_folder,ridtext));
if length(listing) > 1, error('what'); end
if isempty(listing)
    fprintf('\nWarning, no atlas file for rid %d, skipping\n',rid);
    out.enum = [];
    out.ename = {};
    return
end

T = readtable([parcellation_folder,listing.name]);

%% Get the atlas elec labels and clean these
alabels = T.channel;
alabels = decompose_labels(alabels,name);

%% Find the indices of the atlas corresponding to each label
[lia,locb] = ismember(elabels,alabels);

%% Choose the atlas I want
switch which_atlas
    case 'aal'
        anum = T.AAL_region_number;
        aname = T.AAL_label;        
end


%% Assign elecs to atlas
enum = nan(length(elabels),1);
ename = cell(length(elabels),1);
enum(lia) = anum(locb(lia));
ename(lia) = aname(locb(lia));

%% prep output
out.enum = enum;
out.ename = ename;

end