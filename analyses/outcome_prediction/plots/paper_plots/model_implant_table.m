function model_implant_table

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/paper_plots/'];
model_folder = [results_folder,'analysis/outcome/plots/'];

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
% 3rd dimension: mean seeg, std seeg, mean g/s, std g/s, df, tstat, p-value

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
        stat(ia,count,[1,3]) = curr_model.means; % means
        stat(ia,count,[2,4]) = curr_model.sd; % stds
        stat(ia,count,5) = curr_model.stats.df;
        stat(ia,count,6) = curr_model.stats.tstat;
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
        
        fprintf(fid,['Stereo-EEG mean (SD): %1.2f (%1.2f)<br>G/S/D mean (SD):'...
            ' %1.2f (%1.2f)<br><i>t</i>(%d) = %1.2f, %s'],...
            stat(ia,im,1),stat(ia,im,2),...
            stat(ia,im,3),stat(ia,im,4),...
            stat(ia,im,5),stat(ia,im,6),get_p_html(stat(ia,im,7)));
        
        fprintf(fid,'</td>');
        
    end
    
    
    fprintf(fid,'</tr>');
end

%% end table
fprintf(fid,'</table>');
fprintf(fid,'</font>');

fclose(fid);

end