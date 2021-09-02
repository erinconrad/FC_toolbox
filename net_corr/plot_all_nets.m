function plot_all_nets(all)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'corrs/all_net_plots/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

pt_name = all.pt_name;
fn = fieldnames(all.net);
labels = all.labels;
half_indices = 1:2:length(labels);

%% Prep figure
figure
set(gcf,'position',[10 10 1400 500])
tiledlayout(1,length(fn),'tilespacing','tight','padding','tight')

for f = 1:length(fn)
    name = fn{f};
    str = all.net.(name).name;
    net = all.net.(name).data;
    
    
    nexttile
    turn_nans_gray(net)
    
    if f == 1
        yticks(half_indices);
        yticklabels(labels(half_indices))
        
        
    else
        yticklabels([])
        %xticklabels([])
    end
    
    xticks(half_indices);
    xticklabels(labels(half_indices))
    xtickangle(45)
    
    set(gca,'fontsize',15)
    title(str)
    
    if isfield(all.net.(name),'x')
        xlabel(all.net.(name).x);
    end
    
    if isfield(all.net.(name),'y')
        ylabel(all.net.(name).y);
    end
end

print(gcf,[out_folder,pt_name,'_all'],'-dpng')

end