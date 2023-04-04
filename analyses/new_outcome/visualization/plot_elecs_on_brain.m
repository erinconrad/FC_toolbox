function plot_elecs_on_brain(in_name)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = '/data/Human_Data/CNT_iEEG_BIDS/';
inter_folder = [results_folder,'analysis/new_outcome/data/'];
freesurfer_path = '/tools/freesurfer/matlab/';
out_folder = [results_folder,'analysis/new_outcome/plots/elec_locs/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% add freesurfer path
addpath(genpath(freesurfer_path))

%% Load data file
mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;
all_missing = cellfun(@isempty,mt_data.all_spikes(:,1,1));
names = mt_data.all_names;

%% Load Manual validation file
T = readtable('Manual validation.xlsx','Sheet','RIDs');

non_missing = find(~all_missing);
npts = length(non_missing);
for i = 1:npts
    name = names{non_missing(i)};
    if ~isempty(in_name)
        if ~ismember(name,in_name)
            continue
        end
    end

    % find matching row
    r = strcmp(T.name,name);
    assert(sum(r)==1)

    % get rid
    rid = T.RIDs(r);

    % make rid text
    if rid <100
        folder_text = sprintf('sub-RID00%d/',rid);
    else
        folder_text = sprintf('sub-RID0%d/',rid);
    end

    % Load T1 file and get transformation matrix
    t1_file = [data_folder,folder_text,'derivatives/freesurfer/mri/T1.mgz'];
    mri = MRIread(t1_file);
    vox_2_ras = mri.vox2ras;
    tkras = mri.tkrvox2ras;
    transform = @(x) ((vox_2_ras * inv(tkras) * ([x repmat(ones,size(x,1),1)])')');
    first_three_columns = @(x) x(:,1:3);
    all_trans = @(x) first_three_columns(transform(x));

    % load pial files and get vertices and faces and apply transformation
    pial_folder = [data_folder,folder_text,'derivatives/freesurfer/surf/'];
    lobj = SurfStatReadSurf([pial_folder,'lh.pial']);
    lvertices = lobj.coord';
    lfaces = (lobj.tri);

    robj = SurfStatReadSurf([pial_folder,'rh.pial']);
    rvertices = robj.coord';
    rfaces = robj.tri;

    % apply transformation
    rvt = all_trans(rvertices);
    lvt = all_trans(lvertices);

    % Get electrode localizations
    module2 = [data_folder,folder_text,'derivatives/ieeg_recon/module2/']; % module2 folder
    listing = dir([module2,'*.txt']);
    nlist = length(listing);
    found_it = 0;
    for l = 1:nlist
        if contains(listing(l).name,'mm')
            found_it = 1;
            fname = listing(l).name;
            break
        end

    end
    assert(found_it==1)
    T = readtable([module2,fname]);
    locs = [T.Var1 T.Var2 T.Var3];
    allowable_labels = get_allowable_elecs('HUP100');
    mt_symm = find_mt_symmetric_coverage(names,allowable_labels);
    mt = ismember(names,mt_symm);

    % get electrode names
    listing = dir([module2,'*electrode_names.txt']);
    assert(length(listing)==1)
    Tnames = readtable([module2,listing.name],'ReadVariableNames',false);
    names = Tnames.Var1;

    % plot the brain surface
    figure
    rh = trisurf(lfaces,lvt(:,1),lvt(:,2),lvt(:,3));
    hold on
    lh = trisurf(rfaces,rvt(:,1),rvt(:,2),rvt(:,3));
    hold on
    rh.LineStyle = 'none';
    rh.FaceAlpha = 0.1;
    rh.FaceColor = [0.7 0.6 0.6];
    lh.LineStyle = 'none';
    lh.FaceAlpha = 0.1;
    lh.FaceColor = [0.7 0.6 0.6];
    scatter3(locs(~mt,1),locs(~mt,2),locs(~mt,3),'markerfacecolor','k','markeredgecolor','k')
    scatter3(locs(mt,1),locs(mt,2),locs(mt,3),'markerfacecolor','r','markeredgecolor','r')
    
    
    % Name LA, LB, etc.
    end_elecs = {'LA12','LB12','LC12','RA12','RB12','RC12'};
    for e = 1:length(end_elecs)
        curr = end_elecs{e};
        match = strcmp(names,curr);
        if sum(match) ~=0
            curr_loc = locs(match,:);
            if strcmp(curr(1),'L')
                toff = [-10 0 0];
            else
                toff = [10 0 0];
            end
            curr_loc = curr_loc+toff;
            
            if ismember(curr,mt_symm)
                text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','r')
            else
                text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','k')
            end
            
        end
    end
    
    view(-180,-90)%view(-176,4) %view(-182,-5)
    axis off
    print(gcf,[out_folder,name],'-dpng')
    
    close(gcf)


end

end