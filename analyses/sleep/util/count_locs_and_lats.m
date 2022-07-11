function summ = count_locs_and_lats(locs,lats)

possible_locs = {'temporal neocortical','mesial temporal','other cortex'};
possible_lats = {'Left','Right'};
n_possible_locs = length(possible_locs);
n_possible_lats = length(possible_lats);

%% Combine locs and lats to get all possible
possible_loc_lats = cell(n_possible_locs*n_possible_lats,1);
count = 0;
for i = 1:n_possible_locs
    for j = 1:n_possible_lats
        count = count+1;
        possible_loc_lats{count} = [possible_lats{j},' ',possible_locs{i}];
    end
end

npts = length(locs);
loc_lat_count = nan(npts,2);
patient_locs = cell(npts,1);
empty_locs_lats = zeros(npts,2);
alt_lat_nums = zeros(n_possible_lats,2);
for ip = 1:npts
    
    % get their locs and lats
    curr_locs = locs{ip};
    curr_lats = lats{ip};
    
    if all(cellfun(@isempty,curr_locs))
        empty_locs_lats(ip,1) = 1;
        if all(cellfun(@isempty,curr_lats))
            empty_locs_lats(ip,2) = 1;
        end
        continue
    end
    
    % get number of lats
    all_unique_lats = unique(curr_lats);
    matching_lats = cellfun(@(x) any(strcmp(x,possible_lats)),all_unique_lats);
    assert(sum(matching_lats)<=2)
    loc_lat_count(ip,2) = sum(matching_lats);
    
    for il = 1:n_possible_lats
        if any(strcmp(possible_lats{il},all_unique_lats))
            alt_lat_nums(il,1) = alt_lat_nums(il,1)+1;
        end
    end
    
    % combine locs
    combined_locs = cellfun(@(x,y) sprintf('%s %s',x,y), curr_lats,curr_locs,'uniformoutput',false);
    
    % get number of locs
    all_unique_locs = unique(combined_locs);
    matching_locs = cellfun(@(x) any(strcmp(x,possible_loc_lats)),all_unique_locs);
    patient_locs{ip} = all_unique_locs(matching_locs);
    loc_lat_count(ip,1) = sum(matching_locs);
end

%table(loc_lat_count(:,1),loc_lat_count(:,2))
assert(isequal(empty_locs_lats(:,1),empty_locs_lats(:,2)))

alt_lat_nums(:,2) = alt_lat_nums(:,1)/sum(empty_locs_lats(:,1)==0);

%% Get summary on how many patients have coverage of different regions
loc_nums = nan(length(possible_locs),2);
for il = 1:length(possible_locs)
    curr_loc = possible_locs{il};
    loc_nums(il,1) = sum(cellfun(@(x) any(strcmp(x,curr_loc)),locs));
    loc_nums(il,2) = loc_nums(il,1)/sum(empty_locs_lats(:,1)==0);
end

lat_nums = nan(length(possible_lats),2);
for il = 1:length(possible_lats)
    curr_lat = possible_lats{il};
    lat_nums(il,1) = sum(cellfun(@(x) any(strcmp(x,curr_lat)),lats)); 
    lat_nums(il,2) = lat_nums(il,1)/sum(empty_locs_lats(:,1)==0);
end

summ.loc_lat_count = loc_lat_count;
summ.patient_locs = patient_locs;
summ.lat_nums = lat_nums;
summ.loc_nums = loc_nums;
summ.possible_locs = possible_locs;
summ.possible_lats = possible_lats;
summ.empty_locs_lats = empty_locs_lats;

end