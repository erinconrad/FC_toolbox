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
T = readtable([data_folder,'clinical_info/clinical.csv']);

%% Loop over patients in table
for ir = 1:size(T,1)
    %rname = T.ieegportalsubjno{ir};
    rname = T.hupsubjno(ir);
    
    % fix for hup139
    if isnan(rname)
        if T.record_id(ir) == 307
            rname = 139;
        end
    end
    % find the hup***_ string
    
    %{
    
    % see if it contains HUP
    if ~contains(rname,'HUP')
        continue
    end
    
    [si,ei] = regexp(rname,'HUP\d*_');
    
    
    if isempty(si), error('why'); end
    
    rname = rname(si:ei-1);
    %}
    if rname <100
        rname = sprintf('HUP0%d',rname);
    else
        rname = sprintf('HUP%d',rname);
    end
    
    for ip = 1:length(pt)
        name = pt(ip).name;
        if strcmp(name,rname)
            
            % Age 
            pt(ip).clinical.age_implant = T.ageatieegimplant(ir);
            
            % Sex
            sex = T.sex(ir);
            if sex == 1
                pt(ip).clinical.sex = 'Male';
            elseif sex == 2
                pt(ip).clinical.sex = 'Female';
            else
                error('what?');
            end
            
            % Age onset
            pt(ip).clinical.age_onset = T.sz_hist_duration(ir);
            
            % was it stereo?
            %pt(ip).clinical.stereo = logical(T.ieeg_implanttype___4(ir));
            
            % surgery type
            if T.outcome_proctype___1(ir) == 1
                surg = 'Resection';
            elseif T.outcome_proctype___2(ir) == 1
                surg = 'Resection without implant';
            elseif T.outcome_proctype___3(ir) == 1
                surg = 'RNS';
            elseif T.outcome_proctype___4(ir) == 1
                surg = 'Laser ablation';
            elseif T.outcome_proctype___5(ir) == 1
                surg = 'VNS';
            elseif T.outcome_proctype___6(ir) == 1
                surg = 'DBS';
            else
                surg = 'Other';
            end
            pt(ip).clinical.surgery = surg;
            
            % resection hemisphere
            if T.demog_hemi(ir) == 2  
                resection_lat = 'right';
            elseif T.demog_hemi(ir) == 3
                resection_lat = 'left';
            else
                resection_lat = 'na';
            end
            pt(ip).clinical.resection_lat = resection_lat;
            
            % resection type
            if T.demog_resectiontype___1(ir) == 1
                if contains(T.ieeg_resectarea(ir),'temporal lobe') ||...
                        T.demog_resectiontype___2(ir) == 1 || ...
                        T.demog_resectiontype___3(ir) == 1
                    resection_type = 'ATL';
                else
                    resection_type = 'neocortical';
                end
            elseif T.demog_resectiontype___2(ir) == 1
                resection_type = 'ATL';
            elseif T.demog_resectiontype___3(ir) == 1
                if contains(T.ieeg_resectarea(ir),'temporal lobe') || T.demog_resectiontype___2(ir) == 1
                    resection_type = 'ATL';
                else
                    resection_type = 'selective amygdalohippocampectomy';
                end
            elseif T.demog_resectiontype___4(ir) == 1
                resection_type = 'corpus callosotomy';
            elseif T.demog_resectiontype___5(ir) == 1
                resection_type = 'Multiple subpial transection';
            elseif T.demog_resectiontype___6(ir) == 1
                resection_type = 'Radio surgery';
            elseif T.demog_resectiontype___7(ir) == 1
                resection_type = 'Hemispherectomy';
            elseif T.demog_resectiontype___8(ir) == 1
                resection_type = 'therapeutic brain stimulation';
            elseif T.demog_resectiontype___9(ir) == 1
                resection_type = 'other';
            else
                resection_type = 'na';
            end
            pt(ip).clinical.resection_loc = resection_type;
            
            % ablation lateralization
            if T.la_lat_1(ir) == 1
                ablation_lat = 'right';
            elseif T.la_lat_1(ir) == 2
                ablation_lat = 'left';
            else
                ablation_lat = 'na';
            end
            pt(ip).clinical.ablation_lat = ablation_lat;
            
            % Ablation region
            if T.la_region_1___1(ir) == 1
                if T.la_region_1___2(ir) == 1
                    ablation_loc = 'mesial temporal';
                else
                    
                    ablation_loc = 'temporal neocortical';
                end
            elseif T.la_region_1___2(ir) == 1
                ablation_loc = 'mesial temporal';
            elseif T.la_region_1___3(ir) == 1
                ablation_loc = 'frontal';
            elseif T.la_region_1___4(ir) == 1
                ablation_loc = 'insula';
            elseif T.la_region_1___5(ir) == 1
                ablation_loc = 'cingulate';
            elseif T.la_region_1___6(ir) == 1
                ablation_loc = 'parietal';
            elseif T.la_region_1___9(ir) == 1
                ablation_loc = 'other';
            else
                ablation_loc = 'NA';
            end
            pt(ip).clinical.ablation_loc = ablation_loc;
            
            % ILAE outcome
            tilae = [T.demog_ilae1year(ir) T.demog_ilae2years(ir)];
            ilae = cell(2,1);
            for i = 1:length(tilae)
                switch tilae(i)
                    
                    case 1
                        ilae{i} = 'ILAE 1';
                    case 2
                        ilae{i} = 'ILAE 1a';
                    case 3
                        ilae{i} = 'ILAE 2';
                    case 4
                        ilae{i} = 'ILAE 3';
                    case 5
                        ilae{i} = 'ILAE 4';
                    case 6
                        ilae{i} = 'ILAE 5';
                    case 7
                        ilae{i} = 'ILAE 6';
                    otherwise
                        ilae{i} = '';
                end
            end
            ilae_years = [1 2];
            pt(ip).clinical.ilae = ilae;
            pt(ip).clinical.ilae_years = ilae_years;
            
            % engel outcome
            tengel = [T.demog_year(ir) T.demog_years2(ir)];
            engel = cell(2,1);
            engel_years = [1 2];
            for i = 1:length(tengel)
                switch tengel(i)
                    
                    case 1
                        engel{i} = 'IA';
                    case 5
                        engel{i} = 'IB';
                    case 6
                        engel{i} = 'IC';
                    case 7
                        engel{i} = 'ID';
                    case 2
                        engel{i} = 'IIA';
                    case 8
                        engel{i} = 'IIB';
                    case 9
                        engel{i} = 'IIC';
                    case 10
                        engel{i} = 'IID';
                    case 3
                        engel{i} = 'IIIA';
                    case 11
                        engel{i} = 'IIIB';
                    case 4
                        engel{i} = 'IVA';
                    case 12
                        engel{i} = 'IVB';
                    case 13
                        engel{i} = 'IVC';
                    otherwise
                        engel{i} = '';
                        
                end
            end
            pt(ip).clinical.engel = engel;
            pt(ip).clinical.engel_years = engel_years;
            
       
          
            
            
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