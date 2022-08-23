function do_all_figs(doing_from_github)

% close existing figures
close all

%% Parameters
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};


%% Get file locs
locations = fc_toolbox_locs;
plot_folder = locations.paper_plot_folder;

if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
model_folder = locations.paper_plot_folder;

%% Delete the current results file
if exist([plot_folder,'results.html'],'file')~=0
    delete([plot_folder,'results.html']);
end

if exist([plot_folder,'supplemental_results.html'],'file')~=0
    delete([plot_folder,'supplemental_results.html']);
end

%% Load intermediate datasets
% Load results from model
brain_model = load([model_folder,'model_stuff_brainnetome.mat']);
aal_model = load([model_folder,'model_stuff_aal_bernabei.mat']);

% symmetric coverage tests for both atlases
aal_out = load([model_folder,'symm_cov_aal_bernabei.mat']);
aal_out = aal_out.nout;
brain_out = load([model_folder,'symm_cov_brainnetome.mat']);
brain_out = brain_out.nout;

sub_brain = load([model_folder,'nbrain_out.mat']);
sub_brain = sub_brain.nbrain_out;


sub_aal = load([model_folder,'naal_out.mat']);
sub_aal = sub_aal.naal_out;

    %% Supplemental Figure 2 and 3 - electrode subsampling test
subsample_elecs_test(sub_brain,sub_aal,plot_folder)


%% Prep some general results
fid = fopen([plot_folder,'results.html'],'a');
fprintf(fid,['<p>We included all patients who had available electrode localizations (%d patients), although the number '...
    'of patients analyzed varied by analysis, as noted in the results of individual '...
    'analyses. Patients were heterogeneous by age, sex, seizure onset zone localization ',...
    'and lateralization, and implant strategy (Table 1).</p>'],sum(aal_out.pts_with_any_locs));
fclose(fid);

if 1
    


%% Fig 1 - conceptual fig
if doing_from_github == 0
    main_conceptual_figure
end


%% Figure 2 - symmetric coverage test for brainnetome and Supplemental Fig 1 (same but AAL)
% NEED TO THINK ABOUT N
% Prep section in text
symmetric_cov_figure(brain_out,aal_out,plot_folder)



%% Supplemental Fig 4 - coherence (both atlases)
coherence_plots(brain_out,aal_out,plot_folder,freqs)

%% Figure 3 and Figure S5 - confusion matrixes
lat_info(1) = brain_out.lat_info;
lat_info(2) = aal_out.lat_info;
lat_model_figure(lat_info,plot_folder)

%% Table 2 - model comparisons
model_comp_table

%% Figure 4 - null model
null_model_conceptual

%% Figure 5 and Fig S6 - model
model_figure(brain_model,aal_model,plot_folder,doing_from_github)

%% Table S2 - sEEG vs grid/strip/depth implantation
model_implant_table

close all

end