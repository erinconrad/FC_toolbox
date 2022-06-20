function perc_in_gray =count_gray_white(out)

locs = out.circ_out.all_elec_locs;
inc = find(all(out.model_out_gray.excluded==0,2));
ninc = length(inc);
n_gray_white_total = nan(ninc,3);

for ip = 1:ninc
    curr_pt = inc(ip);
    curr_locs = locs{curr_pt};
    
    is_gray = strcmp(curr_locs,'other cortex') | strcmp(curr_locs,'temporal neocortical') | strcmp(curr_locs,'mesial temporal');
    is_white = strcmp(curr_locs,'white matter');
    
    n_gray_white_total(ip,:) = [sum(is_gray),sum(is_white),length(is_gray)];

    
end

perc_in_gray = nanmean(n_gray_white_total(:,1)./n_gray_white_total(:,3));

end