function [ai,signed] = intra_mt_electrode_thing(labels,thing,uni,last_dim)

nelecs = length(labels);
which_elecs = {'A','B','C'};
which_lats = {'L','R'};
maxn = 9;
nmt = length(which_elecs);
intra = nan(nmt,2,last_dim);

% get first two letter in each label
first_two = cellfun(@(x) x(1:2),labels,'uniformoutput',false);

% Get numeric portion
number = cellfun(@(x) (regexp(x,'\d*','Match')),labels,'uniformoutput',false);
number(cellfun(@isempty,number)) = {{'9999'}};
number = cellfun(@(x) str2num(x{1}),number);

for i = 1:nmt
    
    for j = 1:2
    
        curr_elec = [which_lats{j},which_elecs{i}];

        %% Find the contacts matching this electrode
        matching_contacts = strcmp(first_two,curr_elec) & number <= maxn;

        %% Calculate "intra" for these contacts, which depends on what the thing is
        if uni == 1
            curr_intra = nanmean(thing(matching_contacts,:,:),1);
        else
            curr_intra = nanmean(thing(matching_contacts,matching_contacts,:,:),[1 2]);
        end

        %% Fill
        intra(i,j,:) = curr_intra;

    end

end

%% Take AI
ai = abs(intra(:,1,:)-intra(:,2,:))./(intra(:,1,:)+intra(:,2,:));
signed = (intra(:,1,:)-intra(:,2,:))./(intra(:,1,:)+intra(:,2,:));

%% Average across ms
ai = (squeeze(nanmean(ai,1)))';
signed = (squeeze(nanmean(signed,1)))';

end