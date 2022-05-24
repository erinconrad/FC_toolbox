function model_comp_table

%% Get file locs
locations = fc_toolbox_locs;
model_folder = locations.paper_plot_folder;
plot_folder = locations.paper_plot_folder;

for ia = 1:2
    
if ia == 1
    fname = 'Table2.html';
    which_atlas = 'brainnetome';
else
    fname = 'TableS1.html';
    which_atlas = 'aal_bernabei';
end
    
%% Initialize table file
if exist([plot_folder,fname])~=0
    delete([plot_folder,fname]);
end
fid = fopen([plot_folder,fname],'a');

%% Load results from model
model = load([model_folder,'model_stuff_',which_atlas,'.mat']);
model = model.all_out;

%% Get data into table friendly format
mnames = {'Chance','Null model','Null + connectivity','Null + spikes','All'};
nmodels = length(mnames);
stat = nan(nmodels,nmodels);

thing = model.model_p';

% need to fill in 1,.... for all (each one compared to chance - the
% pvalue associated with the AUC)
stat(2,3) = thing(2,3); % density/connectivity
stat(2,4) = thing(2,4); % density/spikes
stat(2,5) = thing(2,5); % density/all
stat(3,4) = thing(3,4); % connectivity/spikes
stat(3,5) = thing(3,5); % connectivity/all
stat(4,5) = thing(4,5); % spikes/all

individ_auc = nan(nmodels,1);
indiv_boot_ci = model.bootci;
for im = 2:nmodels
    individ_auc(im) = mean(model.model_info(im).all_auc);
end


%% Initialize table
fprintf(fid,'<font face="Arial" size="9">');
fprintf(fid,'<table>');

%% First row - variable names
fprintf(fid,'<tr>');

% Make a blank column
fprintf(fid,'<th></th>');

for im = 1:nmodels
    fprintf(fid,'<th>%s</th>',mnames{im});
end
fprintf(fid,'</tr>');

%% FIll in first row with basic significance test of AUC
% Make a column that is the model name
fprintf(fid,'<th>%s</th>',mnames{1});
fprintf(fid,'<th></th>'); % blank one
for im = 2:nmodels
    fprintf(fid,'<td>');
    fprintf(fid,'Mean AUC %1.2f<br>95&#37; CI [%1.2f-%1.2f]',...
        individ_auc(im),indiv_boot_ci(im,1),indiv_boot_ci(im,2));
    fprintf(fid,'</td>');
end

%% Loop over subsequent rows
for im = 2:nmodels
    fprintf(fid,'<tr>');
    
    % Make a column that is the model name
    fprintf(fid,'<th>%s</th>',mnames{im});
    
    for jm = 1:nmodels
        
        
        fprintf(fid,'<td>');
        
        % leave blank if not upper triangle
        if jm <= im
            
        else
            fprintf(fid,'%s',get_p_html(stat(im,jm)));
        end
        
        fprintf(fid,'</td>');
        
    end
    
    fprintf(fid,'</tr>');
end

%% end table
fprintf(fid,'</table>');
fprintf(fid,'</font>');

fclose(fid);

end

end