function rid = simple_hupid_to_rid(id)

%% Get hup id in correct format
if ~strcmp(class(id),'double')
    id = str2num(strrep(id,'HUP',''));
end

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];

%% Load converter
T = readtable([data_folder,'rid_hup.csv']);
hupid = T.hupsubjno;
rids = T.record_id;

r = hupid == id;

if r == 0
    fprintf(sprintf('\nWarning, cannot find hup id %s\n',hup_id))
    rids = nan;
else
    rid = rids(r);
end



end