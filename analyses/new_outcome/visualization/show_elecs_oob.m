function show_elecs_oob

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

figure
for i = 147:length(pt)
    if isempty(pt(i).atropos)
        continue
    end

    % get info
    labels = pt(i).atropos.names;
    locs = pt(i).atropos.xyz;
    atropos = pt(i).atropos.label;
    dkt = pt(i).dkt.label;
    name = pt(i).name;

    % decide if oob
    oob = cellfun(@(x,y) (strcmp(x,'CSF') || strcmp(x,'EmptyLabel')) && strcmp(y,'EmptyLabel'),atropos,dkt);

    % plot
    plot3(locs(:,1),locs(:,2),locs(:,3),'o','color','w');
    hold on
    text(locs(~oob,1),locs(~oob,2),locs(~oob,3),labels(~oob))
    text(locs(oob,1),locs(oob,2),locs(oob,3),labels(oob),'color','r')
    title(name)

    pause
    hold off


end

end