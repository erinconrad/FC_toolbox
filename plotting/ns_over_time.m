function ns_over_time(pc)

%% Parameters
plotm = 1;
spacing = 20;
do_gui = 1;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get network over time
out = net_over_time(pc);
out = reconcile_files(out);


%% Info
ns = out.montage(plotm).ns;
nruns = size(ns,2);
nchs = size(ns,1);
run_center = out.run_center;
file_index = out.file_index;
labels = out.montage(plotm).labels;

%% plot

% raster plot
if do_gui
    figure    
    set(gcf,'position',[10 10 1200 1000])
    turn_nans_gray(ns)
    yticks(1:nchs)
    yticklabels(labels)

    while 1
        try
            [x,~] = ginput;
        catch
            break
        end
        curr_center = run_center(round(x));
        file = file_index(round(x));
        file_name = pc.file(file).name;
        fprintf('\nShowing 15 seconds for %s start time %1.1fs\n',...
            file_name,curr_center-7.5);
        times = [curr_center - 7.5 curr_center + 7.5];
        which_net = 'pc';

        quick_run(file_name,times,which_net,0,3)

    end
end


end