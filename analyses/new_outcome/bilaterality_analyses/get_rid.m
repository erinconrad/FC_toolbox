function rid = get_rid(hup_id)

hup_id = hup_id{1};

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];

%% Load hupid_to_rid csv
T = readtable([data_folder,'atlas/HUPID_to_RID.csv']);
ieegs = T.ieegportalsubjno;
rids = T.record_id;

hup_id_no = strrep(hup_id,'HUP','');
r = contains(ieegs,hup_id);
if r == 0
    temp_hup_id_no = ['0',hup_id_no]; % try adding a leading zero
    temp_hup_id = ['HUP',temp_hup_id_no];
    r = contains(ieegs,temp_hup_id);
end

if r == 0
    if strcmp(hup_id_no(1),'0')
        temp_hup_id_no = hup_id_no(2:end); % try removing leading zero
        temp_hup_id = ['HUP',temp_hup_id_no];
        r = contains(ieegs,temp_hup_id);
    end

end

if r == 0
    fprintf(sprintf('\nWarning, cannot find hup id %s\n',hup_id))
    rid = nan;
    return
end

rid = rids(r);

end