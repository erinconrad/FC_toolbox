function sleep_wake_demographics

%% Parameters

%% Locs
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

%% Get rate in sleep and wake and other variables
rate_sw = bin_out.all_rates;
rate = bin_out.overall_rates;
loc = circ_out.all_locs;
lat = circ_out.all_lats;
sex = circ_out.sex;
age = circ_out.age_implant;
duration = circ_out.duration;

%% Get main demographic info
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
female = strcmp(sex,'Female');
male = strcmp(sex,'Male');
left = strcmp(lat,'left');
right = strcmp(lat,'right');
bilateral = strcmp(lat,'bilateral');

%% Rate in sleep minus rate in wake
rate_diff = (rate_sw(:,2)-rate_sw(:,1));
npts = size(rate_sw,1);

%% Initialize figure
figure

%% Localization
nexttile
plot(1+randn(sum(temporal),1)*0.05,rate_diff(temporal),'o')
hold on
plot(2+randn(sum(extra),1)*0.05,rate_diff(extra),'o')



end

function d = point_to_line(q,p1,p2)

v = (p2-p1)./vecnorm(p2-p1,2,2);
r = (q-p1) - v.*v'.*(q-p1);
d = vecnorm(r,2,2);

end