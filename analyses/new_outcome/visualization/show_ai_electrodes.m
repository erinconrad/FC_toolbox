function show_ai_electrodes(intra,signed,which_elecs,which_lats,name,subplot_path,which_thing,labels)

%% Establish locations for electrodes
nelecs = length(which_elecs);
ncontacts = size(intra,2); % 12
nlats = 2;
locs = repmat(nan(size(intra)),1,1,1,2);
names = cell(size(intra));
mult = 3;
sz = 25;
off = 0.1;
for i = 1:nelecs
    for j = 1:ncontacts
        for k = 1:nlats
            
            names{i,j,k} = [which_lats{k},which_elecs{i},sprintf('%d',j)];

            locs(i,j,k,2) = -i/2;
            
            if k == 1
                locs(i,j,k,1) = -j*mult;
                
            else
                locs(i,j,k,1) = j*mult;
            end

            
            
        end
    end
end

%% Plot the electrodes
% get colormap
cvars_map = make_colormap(512,intra);
cmap = parula(512);

if ~all(isnan(intra),'all')
    figure
    set(gcf,'position',[10 10 1000 600]);
    for i = 1:nelecs
        for j = 1:ncontacts
            for k = 1:nlats
                if ~ismember(names{i,j,k},labels)
                    plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                        'markersize',sz,'markerfacecolor',[1 1 1]);
                else
                    if isnan(intra(i,j,k)) || isnan(cvars_map(i,j,k))
                        plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                            'markersize',sz,'markerfacecolor',[0.7 0.7 0.7]);
                    else
                        
        
                        plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                            'markersize',sz,'markerfacecolor',cmap(cvars_map(i,j,k),:));
                    end

                    hold on
                    text(locs(i,j,k,1),locs(i,j,k,2),sprintf('%1.1f',intra(i,j,k)),...
                        'horizontalalignment','center','fontsize',13);
                    if ismember(names{i,j,k},labels)
                        text(locs(i,j,k,1),locs(i,j,k,2)+off,names{i,j,k},...
                        'horizontalalignment','center','fontsize',13);
                    end

                end
                
                
            end
        end
    end
    ylim([-1.75 -0.25])
    xlim([min(locs(:,:,:,1),[],'all')-2 max(locs(:,:,:,1),[],'all')+2])
    title(sprintf('%s: %s AI %1.1f',name,which_thing,signed),'fontsize',20)
%table(labels,thing)
    %pause
    xticklabels([])
    yticklabels([])
    saveas(gcf,[subplot_path,name,'_',which_thing]);
    close gcf
end



end