function show_locs(p)

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if ischar(p)
    for i = 1:length(pt)
        if strcmp(p,pt(i).name)
            p = i;
            break
        end
    end
end

% Loop over elecs
if isempty(pt(p).elecs)
    fprintf('\nNo electrode info.\n');
    return
end

for i = 1:length(pt(p).elecs)
    locs = pt(p).elecs(i).locs;
    names = pt(p).elecs(i).elec_names;
    
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'w');
    hold on
    text(locs(:,1),locs(:,2),locs(:,3),names,'horizontalalignment','center');
    title(sprintf('Electrode file %d of %d',i,length(pt(p).elecs)))
end

    

end