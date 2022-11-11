function [signed,match] = alt_intra_mt_electrode_thing(labels,thing,uni,last_dim,which_thing)

which_elecs = {'A','B','C'};
which_lats = {'L','R'};

%% Label stuff

%% Misc
% replace '-' with '--'
labels = cellfun(@(x) strrep(x,'-','--'),labels,'uniformoutput',false);

%% Get numeric info
number = cellfun(@(x) (regexp(x,'\d*','Match')),labels,'uniformoutput',false);
% replace empty with really high one
number(cellfun(@isempty,number)) = {{'9999'}};

% Get the first number per label (e.g., LA4--LA5 -> 4)
number = cellfun(@(x) str2num(x{1}),number);

%% Get which electrode
% get the letters
letters = cellfun(@(x) regexp(x,'[a-zA-Z]+','Match'),labels,'uniformoutput',false);
letters(cellfun(@isempty,letters)) = {{'zzzzz'}};
letters = cellfun(@(x) x{1},letters,'uniformoutput',false);
% get first two letter in each label
%first_two = cellfun(@(x) x(1:2),labels,'uniformoutput',false);


maxn = 12; % up to 12 electrodes
nmt = length(which_elecs);


%% Find the matching electrodes
possible_matches = cell(2*nmt*12,1);
count = 0;
for i = 1:nmt
    for k = 1:maxn
        for j = 1:2
            count = count + 1;
            curr = [which_lats{j},which_elecs{i},sprintf('%d',k)];
            possible_matches{count} = curr;
        end
    end
end

match = ismember(possible_matches,labels);
match = possible_matches(match);

switch which_thing{1}

    case {'inter_coh','inter_fc','inter_rl'} % do inter

        % initialize
        inter = nan(nmt,maxn,last_dim);
        
        % loop over electrodes
        for i = 1:nmt
            for k = 1:maxn
                left_elec = strcmp(letters,['L',which_elecs{i}]) & number == k;
                right_elec = strcmp(letters,['R',which_elecs{i}]) & number == k;
    
                if uni == 1
                    % get average rl
                    inter(i,k,:) = nanmean(thing(left_elec|right_elec));

                else
                    % get average inter-electrode connectivity
                    inter(i,k,:) = nanmean(thing(left_elec,right_elec,:),[1 2]);
                end
            end

        end

        % Average across electrodes
        signed = squeeze(nanmean(inter,[1 2]))';

    otherwise % doing intra 

        % initialize
        intra = nan(nmt,maxn,2,last_dim); % n elecs, ncontacts, L and R, bonus

        % which electrodes
        for i = 1:nmt
            % which laterality
            for j = 1:2
                curr_elec = [which_lats{j},which_elecs{i}];
                
                if uni == 1
                    % Can do this on individual contact level
                    for k = 1:maxn
               
                        % Find the contacts matching this electrode
                        matching_contacts = strcmp(letters,curr_elec) & number == k;
                               
                        % Calculate "intra" for these contacts, which depends on what the thing is
                        curr_intra = nanmean(thing(matching_contacts,:,:),1);
                        
                        % Fill
                        intra(i,k,j,:) = curr_intra;
            
                    end
        
                else
        
                    switch which_thing{1}
        
                        case {'near_coh','near_fc'}
                            % measure average FC between one and the one
                            % next to it
                            for k = 1:maxn-1
                                first = strcmp(letters,curr_elec) & number == k;
                                second = strcmp(letters,curr_elec) & number == k+1;
                                curr_intra = nanmean(thing(first,second,:),[1 2]);
                                intra(i,k,j,:) = curr_intra;
        
                            end
        
                        case {'coh','fc'}
                            %% Measure mesial to lateral connectivity
                            
                            mesial_contacts = strcmp(letters,curr_elec) & ismember(number,[1:6]);
                            lateral_contacts = strcmp(letters,curr_elec) & ismember(number,[7:12]);
                
                            % measure connectivity mesial to lateral
                            curr_intra = nanmean(thing(mesial_contacts,lateral_contacts,:),[1 2]);
                
                            
                            % Fill, repeating across all electrodes
                            intra(i,:,j,:) = repmat(curr_intra,1,maxn,1,1);
                    end
                end
        
            end
        
        end
        
        %% Take AI
        %signed = (intra(:,:,1,:)-intra(:,:,2,:))./(intra(:,:,1,:)+intra(:,:,2,:));
        signed = (intra(:,:,1,:)-intra(:,:,2,:))./sqrt((intra(:,:,1,:)).^2+(intra(:,:,2,:)).^2);
        %signed = 0.5-intra(:,:,1,:).*(intra(:,:,2,:))./((intra(:,:,1,:)).^2+(intra(:,:,2,:)).^2);
        
        %% Average across ms
        signed = (squeeze(nanmean(signed,[1 2 3])))';
        
        if last_dim == 1
            signed = nanmean(signed);
        end
end

end