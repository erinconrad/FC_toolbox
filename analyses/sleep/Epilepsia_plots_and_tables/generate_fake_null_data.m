function out = generate_fake_null_data(out)

%{
This function takes the intermediate datasets used to generate the final
statistics and figures, and it makes fake "null" data. The idea is that the
results from this data should be non-significant. Not all analyses have
this fake null data. See below to see which analyses should be affected
when this is called.
%}

%% Randomize rates according to time of day (Figure 2B)
all_tod_rate = out.circ_out.all_tod_rate;
ntimes = size(all_tod_rate,2);
npts = size(all_tod_rate,1);

% shuffle the rates across times
for ip = 1:npts
    % get random perm of times
    p = randperm(ntimes);
    all_tod_rate(ip,:) = all_tod_rate(ip,p);
end

% Refill with fake
out.circ_out.all_tod_rate = all_tod_rate;


%% Randomize sleep/wake spike rates (Figure 2C)
all_rate = out.bin_out.all_rates;
for ip = 1:npts
    p = randperm(2);
    all_rate(ip,:) = all_rate(ip,p);
end
out.bin_out.all_rates = all_rate;

%% Randomize locs (Figure 2D, 4C)
loc = out.circ_out.all_locs;
loc = loc(randperm(npts));
out.circ_out.all_locs = loc;

%% Randomize sleep stage spike rates (figure 2E)
rate_ss = out.seeg_out.rate_ss;
for ip = 1:npts
    p = randperm(5);
    rate_ss(ip,:) = rate_ss(ip,p);
end
out.seeg_out.rate_ss = rate_ss;


%% Randomize spike sequence stuff (Figure 3 A, B)
seq_sw = out.bin_out.seq_sw;
for ip = 1:npts
    curr12 = seq_sw(ip,1:2);
    p = randperm(2);
    seq_sw(ip,1:2) = curr12(p);

    curr34 = seq_sw(ip,3:4);
    p = randperm(2);
    seq_sw(ip,3:4) = curr34(p);
end
out.bin_out.seq_sw = seq_sw;

%% FC Figure 3C
ns_sw = out.bin_out.ns_sw;
for ip = 1:npts
    p = randperm(2);
    ns_sw(ip,:) = ns_sw(ip,p);
end
out.bin_out.ns_sw = ns_sw;

%% Post ictal spike surge (Figure 4A, 4B)
all_pts_spikes_bins = out.sz_out.all_pts_spikes_bins;
ntimes = size(all_pts_spikes_bins,2);
for ip = 1:npts
    % get random perm of times
    p = randperm(ntimes);
    all_pts_spikes_bins(ip,:) = all_pts_spikes_bins(ip,p);
end
out.sz_out.all_pts_spikes_bins = all_pts_spikes_bins;

end