function model_implant_table

%% Get file locs
locations = fc_toolbox_locs;
model_folder = locations.paper_plot_folder;
plot_folder = locations.paper_plot_folder;

fname = 'TableS2.html';
brain_atlas = 'brainnetome';
aal_atlas = 'aal_bernabei';

%% Initialize table file
if exist([plot_folder,fname])~=0
    delete([plot_folder,fname]);
end
fid = fopen([plot_folder,fname],'a');


stat = nan(2,4,7); 
% first dimension which atlas
% second dimension which model
% 3rd dimension: mean seeg, ci95 seeg, mean g/s, ci95 g/s, p-value

atlas_names = {'Brainnetome','AAL'};
mnames = {'Null model','Null + connectivity','Null + spikes','All'};

% loop over atlases
for ia = 1:2
    
    if ia == 1
        which_atlas = brain_atlas;
    else
        which_atlas = aal_atlas;
    end
    
    % Load results from model
    model = load([model_folder,'model_stuff_',which_atlas,'.mat']);
    model = model.all_out;
    
    % Loop over models
    count = 0;
    for im = [2 3 4 5]
        count = count +1;
        curr_model = model.stereo_vs_not.model(im);
        
        % get the stuff
        stat(ia,count,[1,4]) = curr_model.means; % means
        stat(ia,count,2:3) = curr_model.ci95(1,:); % CI-95 seeg
        stat(ia,count,5:6) = curr_model.ci95(2,:); % CI-95 g/s
        stat(ia,count,7) = curr_model.p;
        
    end
    
    
end

%% Initialize table
fprintf(fid,'<font face="Arial" size="9">');
fprintf(fid,'<table>');

%% First row - variable names
fprintf(fid,'<tr>');

% Make a blank column
fprintf(fid,'<th></th>');

for ia = 1:2
    fprintf(fid,'<th>%s</th>',atlas_names{ia});
end
fprintf(fid,'</tr>');

%% Loop over subsequent rows
for im = 1:4
    % next row
    fprintf(fid,'<tr>');
    
    % First column is model name
    fprintf(fid,'<th>%s</th>',mnames{im});
    
    % loop over atlases
    for ia = 1:2
        % next column (cell)
        fprintf(fid,'<td>');
        
        fprintf(fid,['Stereo-EEG mean (95%% CI): %1.2f (%1.2f-%1.2f)<br>G/S/D mean'...
            ' (95%% CI): %1.2f (%1.2f-%1.2f)<br>%s'],...
            stat(ia,im,1),stat(ia,im,2), stat(ia,im,3),...
            stat(ia,im,4),stat(ia,im,5),stat(ia,im,6),...
            get_p_html(stat(ia,im,7)));
        
        fprintf(fid,'</td>');
        
    end
    
    
    fprintf(fid,'</tr>');
end

%% end table
fprintf(fid,'</table>');
fprintf(fid,'</font>');

fclose(fid);

end