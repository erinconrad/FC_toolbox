function [engel1,engel2,ilae1,ilae2,surg_lat] = replace_with_my_outcomes(names,...
    engel1,ilae1,engel2,ilae2,T,soz_locs,surg_lat,surgery)

npts = length(names);
table_names = T.name;
for i = 1:npts
    % find matching row in table
    r = find(strcmp(names{i},table_names));


    %% Check if this is a patient I should have entered manual outcome data for
    if contains(T.name{r},'MP') % only hup
        continue
    end

    if isempty(soz_locs{i})
        continue
    end

    if ~contains(soz_locs{i},'temporal') % only temporal
        continue
    end

    if ~strcmp(surg_lat{i},'right') && ~strcmp(surg_lat{i},'left') % only for unilateral pts
        if ~strcmpi((T.Lat{i}),'right') && ~strcmpi((T.Lat{i}),'left')
            continue
        else
            surg_lat{i} = lower(T.Lat{i});
        end
    end

    if contains(surgery{i},'VNS') || contains(surgery{i},'RNS') || contains(surgery{i},'DBS') % only for resection
        continue
    end

    if strcmp(T.name{r},'HUP178') || strcmp(T.name{r},'HUP122') ...
           || strcmp(T.name{r},'HUP202') % exceptions, did not get surg
        continue
    end

    
    % If made it to here, I should have surgical info
    surg = T.Surg{r};
    if isempty(surg)
         error('why are you missing this outcome?')
    end


    % replace engel and ilae
    engel1(i) = T.x1_yearEngel(r);
    engel2(i) = T.x2_yearEngel(r);

    if isempty(T.x1_yearILAE(r)) || isnan(T.x1_yearILAE(r))
        ilae1{i} = '';
    else
        ilae1{i} = sprintf('ILAE %d',T.x1_yearILAE(r));
    end

    if isempty(T.x2_yearILAE(r)) || isnan(T.x2_yearILAE(r))
        ilae2{i} = '';
    else
        ilae2{i} = sprintf('ILAE %d',T.x2_yearILAE(r));
    end
    

    %T.name(r)
    %T.x1_yearEngel(r)
    %pause

end

end