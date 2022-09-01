function investigate_preimplant

%{
Summary of investigation. Nothing significant by non-parametric unpaired
test, however, I see following trends:

1) Similar number of electrodes implanted regardless of localization
concordance (surprisingly, slightly more electrodes if DISCORDANT
localization)
2) On the other hand, more electrodes (p = 0.09) implanted if discordant
lateralization (makes sense)
3) Higher % of SOZ contacts (p = 0.07) if discordant localization
4) Similar % of SOZ contacts regardless of lateralization concordance
(slightly higher if concordant, but p = 0.38).

And so the trends agree with the model results, although the model had much
more significant results (which could be the result of doing a mixed
effects thing). I am satisfied.
%}

%% File locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
names = out.circ_out.names;

%% pt file as well
pt = load([out_folder1,'pt.mat']);
pt = pt.pt;

%% SOZ
soz = out.bin_out.all_is_soz;

%% NSOZ
nsoz = cellfun(@(x) sum(x==1),soz);

%% Number of electrodes
nelecs = cellfun(@length,soz);

%% Percent soz
perc_soz = cellfun(@(x) sum(x==1)/length(x),soz);

%% Get preimplant data
npts = length(soz);
mri_lesional = cell(npts,1);
concordant_loc = cell(npts,1);
concordant_lat = cell(npts,1);
for ip = 1:npts
    cname = names{ip};
    found_it = 0;
    for jp = 1:length(pt)
        if strcmp(pt(jp).name,cname)
            found_it = 1;
            mri_lesional{ip} = pt(jp).clinical.pre_implant.MRI_lesional;
            concordant_loc{ip} = pt(jp).clinical.pre_implant.concordant_loc;
            concordant_lat{ip} = pt(jp).clinical.pre_implant.concordant_lat;
            break
        end
    end
    assert(found_it == 1)
end

% clean them (turn to 1s and 0s)
mri_lesional = cellfun(@(x) clean_preimplant_designations(x),mri_lesional);
concordant_loc = cellfun(@(x) clean_preimplant_designations(x),concordant_loc);
concordant_lat = cellfun(@(x) clean_preimplant_designations(x),concordant_lat);

%% compare nelecs for those with concordant locs vs not and concordant lats vs not
figure
tiledlayout(1,2)
nexttile
unpaired_plot(nelecs(concordant_loc==1),nelecs(concordant_loc==0),{'Locs concordant','Locs discordant'},'Number electrodes')

nexttile
unpaired_plot(nelecs(concordant_lat==1),nelecs(concordant_lat==0),{'Lats concordant','Lats discordant'},'Number electrodes')

%% compare nsoz for those with concordant locs vs not and concordant lats vs not
figure
tiledlayout(1,2)
nexttile
unpaired_plot(nsoz(concordant_loc==1),nsoz(concordant_loc==0),{'Locs concordant','Locs discordant'},'Number SOZ electrodes')

nexttile
unpaired_plot(nsoz(concordant_lat==1),nsoz(concordant_lat==0),{'Lats concordant','Lats discordant'},'Number SOZ electrodes')

%% compare percent soz for those with concordant locs vs not and concordant lats vs not
figure
tiledlayout(1,2)
nexttile
unpaired_plot(perc_soz(concordant_loc==1),perc_soz(concordant_loc==0),{'Locs concordant','Locs discordant'},'%% SOZ electrodes')

nexttile
unpaired_plot(perc_soz(concordant_lat==1),perc_soz(concordant_lat==0),{'Lats concordant','Lats discordant'},'%% SOZ electrodes')

end