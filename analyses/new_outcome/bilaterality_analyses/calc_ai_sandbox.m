function signed = calc_ai_sandbox(labels,thing,name,labels_plot,uni,last_dim,...
    which_thing,subplot_path,do_plots)

%{
This function takes a list of features for different electrode contacts,
finds the mesial temporal ones, and calculates a single asymmetry index
(for each frequency band) representing the L-R difference.
%}

which_elecs = {'A','B','C'};
which_lats = {'L','R'};
average_level = 'contact'; % side = whole L vs R side; electrode = single electrode; contact = single contact
do_all_intra_elec = 1; % do all contacts for bivariate measures, not just proximal to distal
proximal = 1:6; % proximal contacts are those 1-6
distal = 7:12; % distal contacts are 7-12


%% Label stuff

%% Misc
% replace '-' with '--'
if isempty(labels)
    signed = nan(1,last_dim);
    return
end
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

maxn = 12; % up to 12 contacts per electrode
nmt = length(which_elecs);


%% calculate the AI measurement
switch which_thing{1}

        
    case {'inter_coh','inter_pearson','inter_rl','inter_bp',...
            'inter_plv','inter_xcor','inter_rel_bp','inter_re','lags'} % do inter-electrode (left to right) connectivity (NOT AI)

        % Note lags is up here, because this is a directed measure
        
        % loop over electrodes
        %
        % initialize
        inter = nan(nmt,maxn,last_dim);
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

        
    otherwise % doing AI

        % initialize
        intra = nan(nmt,maxn,2,last_dim); % n elecs, ncontacts, L and R, nfreq
        intra_labels = cell(nmt,maxn,2);

        % which electrodes
        for i = 1:nmt
            
            % which laterality
            for j = 1:2
                curr_elec = [which_lats{j},which_elecs{i}];

                % which of the 12 contacts
                for k = 1:maxn

                    % if this label exists
                    if sum(ismember(labels,sprintf('%s%d',curr_elec,k)))>0
                        intra_labels(i,k,j) = {sprintf('%s%d',curr_elec,k)};
                    end
                end
                
                if isnan(uni) % nelecs

                    % Can do this on individual contact level
                    for k = 1:maxn
                        % Find the contacts matching this electrode
                        matching_contacts = strcmp(letters,curr_elec) & number == k;
                        
                        % intra is the thing
                        intra(i,k,j,:) = nansum(thing(matching_contacts));

                    end

            
                elseif uni == 1
                   
               
                        % Can do this on individual contact level
                    for k = 1:maxn
               
                        % Find the contacts matching this electrode
                        matching_contacts = strcmp(letters,curr_elec) & number == k;

                        % should just have one match?
                        assert(sum(matching_contacts) <= 1)
                               
                        % Calculate "intra" for these contacts, just the
                        % thing for those contacts
                        curr_intra = nanmean(thing(matching_contacts,:,:),1);
                        
                        % Fill
                        intra(i,k,j,:) = curr_intra;
            
                    end
                  
        
                elseif uni == 0 % bivariate measures, will do proximal-to-distal connectivity
        
                    switch which_thing{1}
        
                        case {'near_coh','near_pearson','near_plv'} % (I don't do this)
                            % measure average FC between one and the one
                            % next to it
                            for k = 1:maxn-1
                                first = strcmp(letters,curr_elec) & number == k;
                                second = strcmp(letters,curr_elec) & number == k+1;
                                curr_intra = nanmean(thing(first,second,:),[1 2]);
                                intra(i,k,j,:) = curr_intra;
        
                            end
        
                        case {'coh','pearson','plv','xcor','re'}
                            %% Measure intra-electrode connectivity
                            
                            mesial_contacts = strcmp(letters,curr_elec) & ismember(number,proximal); % first 6 contacts are the proximal contacts
                            lateral_contacts = strcmp(letters,curr_elec) & ismember(number,distal); % last 6 are the distal contacts
                
                            % measure connectivity mesial to lateral
                            if do_all_intra_elec
                                curr_intra = nanmean(thing(strcmp(letters,curr_elec),strcmp(letters,curr_elec),:),[1 2]);
                            else
                                curr_intra = nanmean(thing(mesial_contacts,lateral_contacts,:),[1 2]);
                            end
                      
                            % Fill, repeating across all electrodes
                            intra(i,:,j,:) = repmat(curr_intra,1,maxn,1,1);
                    end
                end
        
            end
        
        end
        
        %% Take AI
        if isnan(uni)
            nleft = sum(intra(:,:,1,:),'all');
            nright = sum(intra(:,:,2,:),'all');
            signed = (nleft-nright)./sqrt(nleft^2+nright^2);
        else
            switch average_level
                case 'contact'
                    % This averages the AI for all L-R contact pairs, which
                    % may overweight the pairs with lower values (e.g., if
                    % very few spikes, the AI may be quite large). May
                    % increase variance?
                    signed = (intra(:,:,1,:)-intra(:,:,2,:))./((intra(:,:,1,:)+intra(:,:,2,:)));
                    signed = (squeeze(nanmean(signed,[1 2 3])))';
                case 'electrode'
                    % This first averages the feature within each
                    % electrode, calculates AI, and then averages across
                    % the three electrodes
                    intra_avg = nanmean(intra,2);
                    signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./((intra_avg(:,:,1,:)+intra_avg(:,:,2,:)));
                    signed = (squeeze(nanmean(signed,[1 2 3])))';
                case 'side'
                    % This averages the feature across the whole side, then
                    % calculates the AI
                    intra_avg = nanmean(intra,[1 2]);
                    signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./((intra_avg(:,:,1,:)+intra_avg(:,:,2,:)));
                    signed = (squeeze(nanmean(signed,[1 2 3])))';
            end
            %}
        end
    
        %% Average across ms  
        if last_dim == 1
            signed = nanmean(signed);
        end
end

%% Plot
if do_plots
    for d = 1:last_dim
        curr_intra = intra(:,:,:,d);
        curr_signed = signed(d);
        curr_thing = sprintf('%s_%d',which_thing{1},d);
        show_ai_electrodes(curr_intra,curr_signed,which_elecs,which_lats,name,subplot_path,curr_thing,labels_plot)
    end
end

%% Error checking
% I would expect that the labels with symmetric coverage should have values
% for intra (unless nan for other reasons)
if 0
    a = permute(intra_labels,[2,3,1]);
    b = permute(intra,[2,3,1]);
    table(a(:,:,1),b(:,:,1))
    table(a(:,:,2),b(:,:,2))
    table(a(:,:,3),b(:,:,3))
    signed
    pause
end

%% Checking nelecs
if 0
    a = permute(intra_labels,[2,3,1]);
    table(a(:,:,1))
    table(a(:,:,2))
    table(a(:,:,3))
    signed
    pause
end

end