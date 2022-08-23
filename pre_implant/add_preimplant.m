function add_preimplant


%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));

%% Get pre-implant data
T = readtable('Manual validation.xlsx','Sheet','Pre-implant data');

%% Loop over rows in the table
for ir = 1:size(T,1)
    
    rname = T.name(ir);
    no_match = 1;
    
    for ip = 1:length(pt)
        name = pt(ip).name;
       
        if strcmp(name,rname)
            no_match = 0;
            pt(ip).clinical.pre_implant.MRI_lesional = T.MRILesional__Y_N_NA__NAMeansNoMRI{ir};
            pt(ip).clinical.pre_implant.concordant_loc = T.ConcordantLocalization_yes_No_Unclear_{ir};
            pt(ip).clinical.pre_implant.concordant_lat = T.ConcordantLateralization_yes_No_Unclear__YesMeansConcordantUnil{ir};
            
        end
    end
    
    if no_match
        error('what')
    end
    
end

%% Look for missing pts
for ip = 1:length(pt)
    if ~isfield(pt(ip).clinical,'pre_implant')
        fprintf('\nWarning, missing pre_implant for %d %s\n',ip,pt(ip).name);
    end
    
end


save([data_folder,'pt.mat'],'pt');

end