function add_demographics

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));

%% Get demographics table
T = readtable([data_folder,'clinical_info/age_sex.csv']);

%% Loop over patients in table
for ir = 1:size(T,1)
    rname = T.ieegportalsubjno{ir};
    
    % see if it contains HUP
    if ~contains(rname,'HUP')
        continue
    end
    
    % find the hup***_ string
    [si,ei] = regexp(rname,'HUP\d*_');
    
    if isempty(si), error('why'); end
    
    rname = rname(si:ei-1);
    
    for ip = 1:length(pt)
        name = pt(ip).name;
        if strcmp(name,rname)
            pt(ip).clinical.age_implant = T.ageatieegimplant(ir);
            sex = T.sex(ir);
            if sex == 1
                pt(ip).clinical.sex = 'Male';
            elseif sex == 2
                pt(ip).clinical.sex = 'Female';
            else
                error('what?');
            end
            pt(ip).clinical.age_onset = T.sz_hist_duration(ir);
            pt(ip).clinical.stereo = T.ieeg_implanttype___4(ir);
        end
    end
    
end

%% Look for missing pts
for ip = 1:length(pt)
    if ~isfield(pt(ip),'clinical')
        fprintf('\nWarning, missing clinical for %d %s\n',ip,pt(ip).name);
    end
    
end

save([data_folder,'pt.mat'],'pt');

end