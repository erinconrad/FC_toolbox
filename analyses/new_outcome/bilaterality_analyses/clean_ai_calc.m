function [AI,left,right] = clean_ai_calc(labels,thing,uni,last_dim)

%{
This function takes a list of features for different electrode contacts,
finds the mesial temporal ones, and calculates a single asymmetry index
(for each frequency band) representing the L-R difference.
%}

which_elecs = {'A','B','C'};
which_lats = {'L','R'};

% this should only happen in the subsampling analysis
if isempty(thing)
    AI = nan(1,last_dim);
    return
end

%% Misc
% replace '-' with '--'
if isempty(labels)
    AI = nan(1,last_dim);
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

%% calculate the feature value for each contact
% initialize
intra = nan(nmt,maxn,2,last_dim); % n elecs, ncontacts, L and R, nfreq

% which electrodes
for i = 1:nmt
    
    % which laterality
    for j = 1:2
        curr_elec = [which_lats{j},which_elecs{i}];

        % which contact
        for k = 1:maxn
   
            % Find the matching contact
            matching_contacts = strcmp(letters,curr_elec) & number == k;

            % should just have one match?
            assert(sum(matching_contacts) <= 1)
   
            % Calculate "intra" for these contacts
            if uni == 1 % if univariate
                curr_intra = nanmean(thing(matching_contacts,:,:),1); %just the thing
            else % if bivariate
                curr_intra = nanmean(thing(matching_contacts,strcmp(letters,curr_elec),:),[1 2]); % average for this contact with all contacts on same electrode
            end
  
            % Fill
            intra(i,k,j,:) = curr_intra;

        end
          
        

    end

end

%% Take AI
% This averages the feature across the whole side, then
% calculates the AI
intra_avg = squeeze(nanmean(intra,[1 2]));
left = intra_avg(1,:);
right = intra_avg(2,:);
AI = (left-right)./(left+right);



end