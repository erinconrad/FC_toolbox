function rids_hup_order

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];


%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
npts = length(pt);

%% Load rid to hup table
T = readtable([data_folder,'rid_hup.csv']);

hup_id = T.hupsubjno;
rid = T.record_id;

%% Prep arrays
out_index = nan(npts,1);
out_hup = cell(npts,1);
out_rid = nan(npts,1);

% Loop over pt
for i = 1:npts
    % get hup name
    name = pt(i).name;

    out_index(i) = i;
    out_hup{i} = name;

    % get numerical portion
    hup_num = str2num(strrep(name,'HUP',''));
    if isempty(hup_num)
        break
    end

    % find the matching index in the lookup table
    row = hup_id == hup_num;
    if sum(row) ~= 1
        fprintf('\nWarning, no rid for %s\n',name);
        continue
    end

    % find the rid
    curr_rid = rid(row);

    % fill up
    
    out_rid(i) = curr_rid;

end

T = table(out_index,out_hup,out_rid)
writetable(T,[data_folder,'ordered_rid_hup.csv'])



end