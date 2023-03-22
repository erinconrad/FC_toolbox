function show_ai_electrodes(intra,signed,which_elecs,which_lats,name,subplot_path,which_thing,labels)
labels = cellfun(@(x) strrep(x,'-Ref',''),labels,'uniformoutput',false);

%% Establish locations for electrodes
nelecs = length(which_elecs);
ncontacts = size(intra,2); % 12
nlats = 2;
locs = repmat(nan(size(intra)),1,1,1,2);
names = cell(size(intra));
mult = 3;
sz = 30;
off = 0.15;
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
    set(gcf,'position',[10 10 1200 300]);
    for i = 1:nelecs
        for j = 1:ncontacts
            for k = 1:nlats
                if ~ismember(names{i,j,k},labels)
                    plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                        'markersize',sz,'markerfacecolor',[1 1 1]);

                    text(locs(i,j,k,1),locs(i,j,k,2)+off,names{i,j,k},...
                        'horizontalalignment','center','fontsize',14);
                else
                    if isnan(intra(i,j,k)) || isnan(cvars_map(i,j,k))
                        plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                            'markersize',sz,'markerfacecolor',[0.7 0.7 0.7]);
                    else
                        
        
                        plot(locs(i,j,k,1),locs(i,j,k,2),'ko','linewidth',2,...
                            'markersize',sz,'markerfacecolor',cmap(cvars_map(i,j,k),:));
                    end

                    hold on
                    text(locs(i,j,k,1),locs(i,j,k,2)-off,sprintf('%1.1f',intra(i,j,k)),...
                        'horizontalalignment','center','fontsize',16);
                    %if ismember(names{i,j,k},labels)
                        text(locs(i,j,k,1),locs(i,j,k,2)+off,names{i,j,k},...
                        'horizontalalignment','center','fontsize',14);
                    %end

                end
                
                
            end
        end
    end
    ylim([-1.75 -0.25])
    xlim([min(locs(:,:,:,1),[],'all')-2 max(locs(:,:,:,1),[],'all')+2])
   % title(sprintf('%s: %s AI %1.1f',name,which_thing,signed),'fontsize',30)
%table(labels,thing)
    
    xticklabels([])
    yticklabels([])
    axis off
    %pause
    saveas(gcf,[subplot_path,name,'_',which_thing]);
    print(gcf,[subplot_path,name,'_',which_thing],'-dpng')
    print(gcf,[subplot_path,name,'_',which_thing],'-depsc')
    close gcf
end



end