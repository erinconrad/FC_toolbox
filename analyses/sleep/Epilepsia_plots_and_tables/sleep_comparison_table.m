function sleep_comparison_table(rate_ss,all_ss,fname)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/sleep/epilepsia/'];



    
%% Initialize table file
if exist([plot_folder,fname])~=0
    delete([plot_folder,fname]);
end
fid = fopen([plot_folder,fname],'a');

%% Get descriptive stats
nss = length(all_ss);
ss_median= nan(nss,1);
ss_iqr = nan(nss,2);

for is = 1:nss
    ss_median(is) = nanmedian(rate_ss(:,is));
    ss_iqr(is,:) = [prctile(rate_ss(:,is),25) prctile(rate_ss(:,is),75)];
end

%% Do post-hoc sign rank tests
p_ss = nan(nss,nss);
for is = 1:nss
    for js = 1:is-1
        p = signrank(rate_ss(:,is),rate_ss(:,js));
        p_ss(is,js) = p;
    end
end
p_ss = p_ss';


%% Initialize table
fprintf(fid,'<font face="Arial" size="9">');
fprintf(fid,'<table>');

%% First row - variable names
fprintf(fid,'<tr>');

% Make a blank column
fprintf(fid,'<th></th>');

% Make a blank column
%fprintf(fid,'<th></th>');

for im = 1:nss
    fprintf(fid,'<th>%s</th>',all_ss{im});
end
fprintf(fid,'</tr>');

%% FIll in first row with basic significance test of AUC
% Make a column that is the model name
fprintf(fid,'<th>%s</th>','Descriptive statistics');
%fprintf(fid,'<th></th>'); % blank one
for im = 1:nss
    fprintf(fid,'<td>');
    fprintf(fid,'Median: %1.2f, IQR %1.2f-%1.2f',...
        ss_median(im),ss_iqr(im,1),ss_iqr(im,2));
    fprintf(fid,'</td>');
end

%% Loop over subsequent rows
for im = 1:nss
    fprintf(fid,'<tr>');
    
    % Make a column that is the model name
    fprintf(fid,'<th>%s</th>',all_ss{im});

    %fprintf(fid,'<th></th>'); % blank one
    
    for jm = 1:nss
        
        
        fprintf(fid,'<td>');
        
        % leave blank if not upper triangle
        if jm <= im
            
        else
            fprintf(fid,'%s',get_p_html(p_ss(im,jm)));
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