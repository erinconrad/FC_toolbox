function txf_rl = transform_rl(rl,rate,rate_min)

if all(isnan(rate)), txf_rl = nan(size(rl)); return; end

% set rl for those without spikes to inf
rl(rate<rate_min) = inf;

% rank electrodes by rl
[~,I] = sort(rl);
r = (1:length(rl))';
r(I) = r;

% transform ranking to a 0-1 scale (thus normalizing for number of
% electrodes)
txf_rl = map_numbers_onto_range(r,[0 1]);

if 0
    table(rl,r,txf_rl)
end

end