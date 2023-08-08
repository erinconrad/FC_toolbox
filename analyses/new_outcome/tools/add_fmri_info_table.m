function T = add_fmri_info_table(T,fT,rid_table,csv_path)

%% Define temporal ROIs
temporal_hippo_amygdala_left = [108 110 112 114 116 118  74  78  86 212 214 216];
temporal_hippo_amygdala_right = [109 111 113 115 117 119  75  79  87 213 215 217];

% Add one to the indices of the regions (python to matlab)
temporal_hippo_amygdala_left = temporal_hippo_amygdala_left + 1;
temporal_hippo_amygdala_right = temporal_hippo_amygdala_right + 1;

%% Add new variable to T
T = addvars(T,nan(size(T,1),1),'NewVariableNames','fmri_AI');
T = addvars(T,nan(size(T,1),1),'After','names','NewVariableNames','rid');

% Loop over patients in T
for it = 1:size(T,1)
    % get hup ID string
    name = T.names{it};

    % get numeric portion
    hupid = regexp(name,'\d*','Match');
    if isempty(hupid), continue; end
    hupid = str2num(hupid{1});

    % get the RID
    rid = convert_rid_hupid(hupid,'hupid',rid_table);
    

    % skip if you can't find it
    if isempty(rid) || isnan(rid)
        continue
    end

    % add it to table
    T.rid(it) = rid;

    % put it into the format the fmri table expects
    rid_string = sprintf('sub-RID0%d',rid);

    % find the correct row of the fmri table
    fmri_row = strcmp(fT.Subject,rid_string);
    assert(sum(fmri_row) <= 1)

    % skip if you can't find it
    if sum(fmri_row) == 0, continue; end

    % Load the csv containing the patient-specific fcon
    fcon_T = readtable([csv_path,rid_string,'.csv']);
    fcon = table2array(fcon_T);

    % Get the left and right strength
    left_str = mean(sum(abs(fcon(temporal_hippo_amygdala_left,temporal_hippo_amygdala_left)),2),1);
    right_str = mean(sum(abs(fcon(temporal_hippo_amygdala_right,temporal_hippo_amygdala_right)),2),1);

    % Define AI
    fmri_AI = (left_str-right_str)./(left_str+right_str);

    % Add this as a column to the T table
    T.fmri_AI(it) = fmri_AI;


end

if 0
    oT = table(T.names,T.rid,T.fmri_AI,'VariableNames',{'HUPID','RID','AI'});
    oT(strcmp(T.soz_locs,'temporal') & contains(T.names,'HUP'),:)
end

% double checking some things
any_fmri_ieeg = strcmp(fT.IEEG,'IEEG');
rid_text_fmri_eeg = fT.Subject(any_fmri_ieeg);
rid_fmri_eeg = cellfun(@(x) str2num(x(end-2:end)),rid_text_fmri_eeg);
hupid_fmri_eeg = arrayfun(@(x) convert_rid_hupid(x,'rid',rid_table),rid_fmri_eeg);

ooT = table(rid_text_fmri_eeg,hupid_fmri_eeg);

end