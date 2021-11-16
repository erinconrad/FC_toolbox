function output_all_sz_info

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];
out_folder = [locations.main_folder,'data/sz_times/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

T = array2table(zeros(0,7),'VariableNames',...
    {'Patient','IEEGID','IEEGname','start','end','source','notes'});

for p = 1:length(pt)
    
    for f = 1:length(pt(p).ieeg.file)
        if isfield(pt(p).ieeg.file(f),'sz_times')
            for s = 1:size(pt(p).ieeg.file(f).sz_times,1)
                name = pt(p).name;
                fid = f;
                ieegname = pt(p).ieeg.file(f).name;
                sz_start = pt(p).ieeg.file(f).sz_times(s,1);
                sz_end = pt(p).ieeg.file(f).sz_times(s,2);
                source = pt(p).ieeg.file(f).sz_time_source;
                
                if strcmp(source,'portal times, worst source')
                    note = 'left pad 60 s, right pad 10 minutes';
                elseif strcmp(source,'Erin Excel')
                    note = 'looked for seizure-y annotations and portal times and verified onset plus or minus several seconds';
                else
                    note = '';
                end
                
                % add row
                T = [T;{name,fid,ieegname,sz_start,sz_end,source,note}];
                
            end
        end
        
    end
    
    T = [T;{'',nan,'',nan,nan,'',''}];
    
end

writetable(T,[out_folder,'all_szs_table.csv']);

end