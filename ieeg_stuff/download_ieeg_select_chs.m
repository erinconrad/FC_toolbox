function data = download_ieeg_select_chs(fname, login_name, pwfile, ...
    run_times,chs_to_download,name)

attempt = 1;

% Wrap data pulling attempts in a while loop
while 1
    
    try
        session = IEEGSession(fname, login_name, pwfile);
        channelLabels = session.data.channelLabels;
        cl_labels = decompose_labels(channelLabels,name);

        % get fs
        data.fs = session.data.sampleRate;
        
        % find the matching chs
        ch_idx = ismember(cl_labels,chs_to_download);
        ch_idx = find(ch_idx);

        % Convert times to indices
        run_idx = round(run_times(1)*data.fs):round(run_times(2)*data.fs);

        if ~isempty(run_idx)
            % Break the number of channels in half to avoid wacky server errors
            values = session.data.getvalues(run_idx,ch_idx);
            
        else
            values = [];
        end

        data.values = values;

        % get file name
        data.file_name = session.data.snapName;

        % Get ch labels
        data.chLabels = channelLabels(ch_idx);
       
        % break out of while loop
        break
        
    % If server error, try again (this is because there are frequent random
    % server errors).
    catch ME
        if contains(ME.message,'503') || contains(ME.message,'504') || ...
                contains(ME.message,'502') || contains(ME.message,'500')
            attempt = attempt + 1;
            fprintf('Failed to retrieve ieeg.org data, trying again (attempt %d)\n',attempt); 
        else
            ME
            error('Non-server error');
            
        end
        
    end
end

%% Delete session
session.delete;
clearvars -except data

end