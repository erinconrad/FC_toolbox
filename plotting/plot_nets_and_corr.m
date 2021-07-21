function plot_nets_and_corr(out,im)

do_r2 = 1;

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
    dataw = out.montage(im).net(in).data;

    % Expand
    data = wrap_or_unwrap_adjacency_fc_toolbox(dataw);
    

    if do_r2
        data = (data).^2;
    end

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
    if do_r2
        [r,p] = corr((out.montage(im).net(1).data).^2,...
            (out.montage(im).net(2).data).^2,'rows','pairwise');
    else
        [r,p] = corr(out.montage(im).net(1).data,out.montage(im).net(2).data,'rows','pairwise');
    end
    
    if in == 1
        ylabel(sprintf('%s\nr = %1.2f, p = %1.3f',montage,r,p),'fontsize',15);
    end


    if im == 2
        xlabel(network,'fontsize',15)
    end

end



  

end