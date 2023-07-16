function [signed,left,right] = calc_ai_ns(labels,thing,name,labels_plot,atropos,dkt,uni,last_dim,...
    which_thing,subplot_path,do_plots)

%{
This function takes a list of features for different electrode contacts,
finds the mesial temporal ones, and calculates a single asymmetry index
(for each frequency band) representing the L-R difference.
%}
attempt_rm_oob = 0;
which_elecs = {'A','B','C'};
which_lats = {'L','R'};
average_level = 'side'; % side = whole L vs R side; electrode = single electrode; contact = single contact

% this should only happen in the subsampling analysis
if isempty(thing)
    signed = nan(1,last_dim);
    return
end

%% Misc
% replace '-' with '--'
if isempty(labels)
    signed = nan(1,last_dim);
    return
end
labels = cellfun(@(x) strrep(x,'-','--'),labels,'uniformoutput',false);

%% Rule to remove elecs not in brain
if isempty(atropos)
    out_of_brain = zeros(length(atropos),1);
else
    out_of_brain = cellfun(@(x,y) (strcmp(x,'CSF') || strcmp(x,'EmptyLabel')) && strcmp(y,'EmptyLabel'),atropos,dkt);
end



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
letter_no_side = cellfun(@(x) x(2),letters,'uniformoutput',false);

maxn = 12; % up to 12 contacts per electrode
nmt = length(which_elecs);

%% Remove out of brain and contralateral
letter_no_side_number = cellfun(@(x,y) sprintf('%s%d',x,y),letter_no_side,num2cell(number),'uniformoutput',false);

% find the symmetric pair for each one
opp_partner = nan(length(letter_no_side_number),1);
for i = 1:length(letter_no_side_number)
    poss_opp_partners = find(strcmp(letter_no_side_number,letter_no_side_number(i)));
    poss_opp_partners(poss_opp_partners==i) = [];
    if length(poss_opp_partners) == 0
        error('what')
    end
    opp_partner(i) = poss_opp_partners(1);
end


if attempt_rm_oob
    % Define ones to remove
    to_rm_idx = [find(out_of_brain),opp_partner(out_of_brain)];
    to_rm = zeros(length(letter_no_side_number),1);
    to_rm(to_rm_idx) = 1;
    to_rm = logical(to_rm);
    
    if sum(to_rm) >0 % I would not expect any of these to remain.
        error('what')
    end

    if 0
        table(atropos,dkt,out_of_brain,letter_no_side_number,to_rm)
    end
    
    % Make the thing nans
    thing(to_rm) = nan;
end

%% calculate the AI measurement
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
        

        for k = 1:maxn
   
            % Find the contacts matching this electrode
            matching_contacts = strcmp(letters,curr_elec) & number == k;

            % should just have one match?
            assert(sum(matching_contacts) <= 1)

            if contains(which_thing,'inter') % some inter-electrode measure
                if uni == 1
                    curr_intra = nanmean(thing(matching_contacts,number <= maxn,:),1); %just the thing
                else
                    inter_match = strcmp(letter_no_side,which_elecs{i} & number <= maxn);
                    curr_intra = nanmean(thing(matching_contacts,inter_match,:),[1 2]); %just the thing
                end


            else
                   
                % Calculate "intra" for these contacts
                if uni == 1 % if univariate
                    curr_intra = nanmean(thing(matching_contacts,:,:),1); %just the thing
                else % if bivariate
                    curr_intra = nanmean(thing(matching_contacts,strcmp(letters,curr_elec),:),[1 2]); % average for this contact with all contacts on same electrode
                end

            end
            
            
            % Fill
            intra(i,k,j,:) = curr_intra;

        end
          
        

    end

end

%% Take AI
if contains(which_thing,'inter')
    signed = squeeze(nanmean(intra,[1 2 3]))';
else
    switch average_level
        case 'contact'
            % This averages the AI for all L-R contact pairs, which
            % may overweight the pairs with lower values (e.g., if
            % very few spikes, the AI may be quite large). May
            % increase variance?
            %signed = (intra(:,:,1,:)-intra(:,:,2,:))./sqrt((intra(:,:,1,:).^2+intra(:,:,2,:).^2));
            signed = (intra(:,:,1,:)-intra(:,:,2,:))./((intra(:,:,1,:)+intra(:,:,2,:)));
            signed = (squeeze(nanmean(signed,[1 2 3])))';
        case 'electrode'
            % This first averages the feature within each
            % electrode, calculates AI, and then averages across
            % the three electrodes
            intra_avg = nanmean(intra,2);
            %error('what')
            %signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./sqrt((intra_avg(:,:,1,:).^2+intra_avg(:,:,2,:).^2));
            signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./(intra_avg(:,:,1,:)+intra_avg(:,:,2,:));
            signed = (squeeze(nanmean(signed,[1 2 3])))';
        case 'side'
            % This averages the feature across the whole side, then
            % calculates the AI
            intra_avg = nanmean(intra,[1 2]);
            %error('what')
            %signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./(sqrt(intra_avg(:,:,1,:).^2+intra_avg(:,:,2,:).^2));
            signed = (intra_avg(:,:,1,:)-intra_avg(:,:,2,:))./(intra_avg(:,:,1,:)+intra_avg(:,:,2,:));
            signed = (squeeze(nanmean(signed,[1 2 3])))';
            left = nanmean(intra_avg(:,:,1,:));
            right = nanmean(intra_avg(:,:,2,:));
    end
end


%% Average across ms  
if last_dim == 1
    signed = nanmean(signed);
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