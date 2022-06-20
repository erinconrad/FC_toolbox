function missing_stuff(pt,out)

npts = length(pt);
not_stereo = nan(npts,1);
names = cell(npts,1);
missing_locs = nan(npts,1);

for ip = 1:npts
    names{ip} = pt(ip).name;
    not_stereo(ip) = ~pt(ip).clinical.stereo;
    missing_locs(ip) = isempty(pt(ip).elecs);
end

fprintf('\nPatient structure:\n');
table(names,not_stereo,missing_locs)
fprintf('\n%d grid/strip patients, %d with locations.\n',sum(not_stereo),sum(not_stereo&~missing_locs));


npts_alt = length(out.all_names);
not_stereo_alt = nan(npts_alt,1);
names_alt = cell(npts_alt,1);
missing_locs_alt = nan(npts_alt,1);
for ip = 1:npts_alt
    names_alt{ip} = out.all_names{ip};
    not_stereo_alt(ip) = ~out.all_stereo(ip);
    missing_locs_alt(ip) = all(isnan(out.all_locs{ip}),'all');
end

fprintf('\nOut structure:\n');
table(names_alt,not_stereo_alt,missing_locs_alt)
fprintf('\n%d grid/strip patients, %d with locations.\n',sum(not_stereo_alt),sum(not_stereo_alt&~missing_locs_alt));
end