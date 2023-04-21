function plot_elecs_on_brain(overwrite,in_name)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = '/data/Human_Data/CNT_iEEG_BIDS/';
inter_folder = [results_folder,'analysis/new_outcome/data/'];
freesurfer_path = '/tools/freesurfer/matlab/';
out_folder = [results_folder,'analysis/new_outcome/plots/elec_locs/'];
other_data_folder = [locations.main_folder,'data/'];
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

%% Load pt file
pt = load([other_data_folder,'pt.mat']);
pt = pt.pt;

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

    if exist([out_folder,name,'.fig'],'file') ~=0
        if overwrite == 0
            fprintf('\nSkipping %s\n',name);
            continue
        end
    end
    fprintf('\nDoing %s\n',name);
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
    if ~exist(t1_file,'file')
        fprintf('\nNo T1 for %s, skipping\n',name);
        continue
    end
    mri = MRIread(t1_file);
    vox_2_ras = mri.vox2ras;
    tkras = mri.tkrvox2ras;
    transform = @(x) ((vox_2_ras * inv(tkras) * ([x repmat(ones,size(x,1),1)])')');
    first_three_columns = @(x) x(:,1:3);
    all_trans = @(x) first_three_columns(transform(x));

    % load pial files and get vertices and faces and apply transformation
    pial_folder = [data_folder,folder_text,'derivatives/freesurfer/surf/'];

    if ~exist(pial_folder,'dir')
        fprintf('\nNo pial dir for %s, skipping\n',name);
        continue
    end

    if ~exist([pial_folder,'lh.pial'],'file')
        fprintf('\nNo pial file for %s, skipping\n',name);
        continue
    end

    lobj = SurfStatReadSurf([pial_folder,'lh.pial']);
    lvertices = lobj.coord';
    lfaces = (lobj.tri);

    robj = SurfStatReadSurf([pial_folder,'rh.pial']);
    rvertices = robj.coord';
    rfaces = robj.tri;

    % apply transformation
    rvt = all_trans(rvertices);
    lvt = all_trans(lvertices);

    % get electrode names
    module2 = [data_folder,folder_text,'derivatives/ieeg_recon/module2/']; % module2 folder
    listing = dir([module2,'*electrode_names.txt']);
    if length(listing)==0
        fprintf('\nNo elec names for %s, skipping\n',name);
        continue
    end
    Tnames = readtable([module2,listing.name],'ReadVariableNames',false);
    enames = Tnames.Var1;

    % Get electrode localizations
    
    listing = dir([module2,'*.txt']);
    nlist = length(listing);
    found_it = 0;
    for l = 1:nlist
        if contains(listing(l).name,'mm_electrodes')
            found_it = 1;
            fname = listing(l).name;
            break
        end

    end
    
    if found_it == 0
        fprintf('\nNo elec locs for %s, skipping\n',name);
        continue
    end
    eT = readtable([module2,fname]);
    locs = [eT.Var1 eT.Var2 eT.Var3];
    allowable_labels = get_allowable_elecs('HUP100');
    mt_symm = find_mt_symmetric_coverage(enames,allowable_labels);
    mt = ismember(enames,mt_symm);

    
    % get the atlas parcellations
    found_it = 0;
    for p = 1:length(pt)
        if strcmp(name,pt(p).name)
            found_it = 1;
            break
        end
    end
    assert(found_it == 1)

    if isempty(pt(p).atropos)
        atropos_labels = [];
        atropos_names = [];
        dkt_labels = [];
        dkt_names = [];
    else
        atropos_labels = pt(p).atropos.label;
        atropos_names = pt(p).atropos.names;
        dkt_labels = pt(p).dkt.label;
        dkt_names = pt(p).dkt.names;
    end
    
    assert(isequal(dkt_names,atropos_names))

    [Lia,Locb] = ismember(enames,atropos_names);
    assert(isequal(atropos_names(Locb(Lia)),enames(Lia)))
    atropos_labels(Lia) = atropos_labels(Locb(Lia)); atropos_labels(~Lia) = {''};
    dkt_labels(Lia) = dkt_labels(Locb(Lia)); dkt_labels(~Lia) = {''};

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
    scatter3(locs(~mt,1),locs(~mt,2),locs(~mt,3),'markerfacecolor','w','markeredgecolor','k')
    scatter3(locs(mt,1),locs(mt,2),locs(mt,3),'markerfacecolor','w','markeredgecolor','r')
    text(locs(:,1),locs(:,2),locs(:,3),...
        cellfun(@(x,y,z) sprintf('%s, %s, %s',x,y,z), enames, atropos_labels, dkt_labels,'uniformoutput',false),...
        'HorizontalAlignment','center','fontsize',8)
    
    
    %{
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
    %}
    
    view(-180,-90)%view(-176,4) %view(-182,-5)
    axis off
    print(gcf,[out_folder,name],'-dpng')
    savefig(gcf,[out_folder,name,'.fig'])
    
    close(gcf)


end

end