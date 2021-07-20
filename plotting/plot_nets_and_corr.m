function plot_nets_and_corr(out,im)

%% initialize figure
figure
set(gcf,'position',[20 20 1300 700])
tiledlayout(1,2)

% Loop over montages

montage = out.montage(im).name;
labels = out.montage(im).labels;
is_run = out.montage(im).is_run;

% remove labels not in run
labels = labels(is_run);

for in = 1:2
    network = out.montage(im).net(in).name;
    data = out.montage(im).net(in).data;

    % Expand
    data = wrap_or_unwrap_adjacency_fc_toolbox(data);


    % remove things not in run
    data = data(is_run,is_run);


    % show network
    nexttile
    turn_nans_gray(data)
    xticks(1:length(labels))
    yticks(1:length(labels))
    xticklabels(labels)
    yticklabels(labels)
    set(gca,'fontsize',10)

    % calculate correlation between wrapped matrices
    [r,p] = corr(out.montage(im).net(1).data,out.montage(im).net(2).data,'rows','pairwise');
    
    if in == 1
        ylabel(sprintf('%s\nr = %1.2f, p = %1.3f',montage,r,p),'fontsize',15);
    end


    if im == 2
        xlabel(network,'fontsize',15)
    end

end



  

end