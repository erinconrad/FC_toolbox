function outT = five_sense_calculator(impute_missing_data,which_pts)

if exist('impute_missing_data','var') == 0
    impute_missing_data = 0;
end

%% Locations
%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load in manual validation preclinical table
T = readtable('Manual validation.xlsx','Sheet','Pre-implant data');

%% Prep variables
names = T.name;
npts = length(names);
lesion = cell(npts,1);
interictal = cell(npts,1);
neuropsych = cell(npts,1);
semiology = cell(npts,1);
ictal = cell(npts,1);
prob = nan(npts,1);
final_focal = cell(npts,1);
lesion_coefficient = nan(npts,1);
interictal_coefficient = nan(npts,1);
neuropsych_coefficient = nan(npts,1);
semiology_coefficient = nan(npts,1);
ictal_coefficient = nan(npts,1);
intercept = nan(npts,1);

%% Map manual validation entries onto 5-SENSE categories
for i = 1:npts
    name = names{i};
    row = strcmp(T.name,name);
    assert(sum(row)==1)

    % Lesion extent - per paper, this only looks at MRI lesion and can be
    % no lesion, focal lesion, or all other lesions
    mri_lesional_yn = T.MRILesional__Y_N_NA__NAMeansNoMRI{row};
    mri_lat = T.MRILesionalLaterality_left_Right_Bilateral_Broad_NA__NAMeansNon{row};
    mri_loc = T.MRILesionalLocalization_temporal_Frontal_Other_Multifocal_Broad{row};

    if isempty(mri_lesional_yn)
    elseif strcmpi(mri_lesional_yn,'N')
        lesion{i} = 'no lesion';
    else % there is a lesion
        if strcmpi(mri_lat,'bilateral') || strcmpi(mri_lat,'multifocal')
            lesion{i} = 'all other lesions';
        elseif strcmpi(mri_loc,'multifocal') || strcmpi(mri_loc,'broad')
            lesion{i} = 'all other lesions';
        else % focal
            lesion{i} = 'focal';
        end
        
    end

    % interictal EEG extent
    spike_lat = T.ScalpEEGSpikeLaterality_left_Right_Bilateral_Broad_NA__NAMeansN{row};
    spike_loc = T.ScalpEEGSpikeLocalization_temporal_Frontal_Other_Multifocal_Bro{row};

    if isempty(spike_lat)
    elseif strcmpi(spike_lat,'NA')
        interictal{i} = 'no IEDs';
    elseif strcmpi(spike_lat,'bilateral') || strcmpi(spike_lat,'broad') || strcmpi(spike_lat,'unclear')
        interictal{i} = 'bilateral';
    else
        interictal{i} = 'all others';
    end

    % neuropsych deficit
    neuropsych_table = T.neuropsychLanguageDysfunctionLateralization_left_Right_Bilatera{row};

    if isempty(neuropsych_table)
    elseif strcmpi(neuropsych_table,'NA') % means no neuropsych done
    elseif strcmpi(neuropsych_table,'unclear')
        neuropsych{i} = 'no deficit';
    elseif strcmpi(neuropsych_table,'bilateral')
        neuropsych{i} = 'not localizing';
    else
        neuropsych{i} = 'localizing';
    end

    % semiology - per paper, strong vs weak or not localizing
    sem_lat = T.SemiologyLateralization_left_Right_Bilateral_Unclear_NA__NAMean{row};
    sem_loc = T.SemiologyLocalization_temporal_Frontal_Other_Multifocal_Broad_u{row};

    if isempty(sem_lat)
    elseif strcmpi(sem_lat,'NA')
    elseif (strcmpi(sem_lat,'unclear') || strcmpi(sem_lat,'bilateral')) ...
            && (strcmpi(sem_loc,'unclear') || strcmpi(sem_loc,'broad') || strcmpi(sem_loc,'multifocal'))
        % if both lat and loc unclear
        semiology{i} = 'weak or not localizing';
    else
        semiology{i} = 'strong';
    end

    % ictal EEG - per paper, focal/lobar, bilateral synchronous,
    % none/multilobar/diffuse - we've never described bilateral
    % synchronous.
    sz_lat = T.ScalpEEGSeizuresLaterality_left_Right_Bilateral_BroadOrUnclear_{row};
    sz_loc = T.ScalpEEGSeizureLocalization_temporal_Frontal_Other_Multifocal_B{row};

    if isempty(sz_lat)
    elseif strcmpi(sz_lat,'NA') || strcmpi(sz_loc,'NA')
        ictal{i} = 'none/multilobar/diffuse';
    elseif strcmpi(sz_lat,'unclear') || strcmpi(sz_lat,'bilateral') || ...
            strcmpi(sz_lat,'broad') || strcmpi(sz_loc,'multifocal') || ...
            strcmpi(sz_loc,'broad') || strcmpi(sz_loc,'unclear')
        ictal{i} = 'none/multilobar/diffuse';
    else
        ictal{i} = 'focal';
    end

    % Final post-implant determinationf of focality
    final_focal_table = T.Post_implantSeizureOnsetDetermination_unifocal_Bifocal_multifoc{row};
    if isempty(final_focal_table) || strcmpi(final_focal_table,'NA')
    elseif contains(final_focal_table,'unifocal','ignorecase',true)
        final_focal{i} = 'unifocal';
    elseif strcmpi(final_focal_table,'bifocal') || strcmpi(final_focal_table,'broad') || ...
            strcmpi(final_focal_table,'multifocal') || strcmpi(final_focal_table,'missed')
        final_focal{i} = 'other';
    end


    %% Now, apply the model coefficients to get the 5-SENSE score
    intercept(i) = -0.3135; % from paper appendix

    % lesion coefficient
    if isempty(lesion{i})
        lesion_coefficient(i) = nan;
    elseif strcmp(lesion{i},'no lesion')
        lesion_coefficient(i) = -2.2626;
    elseif strcmpi(lesion{i},'focal')
        lesion_coefficient(i) = 0;
    elseif strcmp(lesion{i},'all other lesions')
        lesion_coefficient(i) = -2.1494;
    else
        error('why')
    end

    % interictal EEG coefficient
    if isempty(interictal{i})
        interictal_coefficient(i) = nan;
    elseif strcmp(interictal{i},'no IEDs')
        interictal_coefficient(i) = 1.8056;
    elseif strcmp(interictal{i},'bilateral')
        interictal_coefficient(i) = 0;
    elseif strcmp(interictal{i},'all others')
        interictal_coefficient(i) = 1.1807;
    else
        error('why')
    end

    % neuropsych coefficient
    if isempty(neuropsych{i})
        neuropsych_coefficient(i) = nan;
    elseif strcmp(neuropsych{i},'no deficit')
        neuropsych_coefficient(i) = 1.1550;
    elseif strcmp(neuropsych{i},'not localizing')
        neuropsych_coefficient(i) = -0.2554;
    elseif strcmp(neuropsych{i},'localizing')
        neuropsych_coefficient(i) = 0;
    else
        error('why')
    end

    % semiology coefficient
    if isempty(semiology{i})
        semiology_coefficient(i) = nan;
    elseif strcmp(semiology{i},'weak or not localizing')
        semiology_coefficient(i) = 0;
    elseif strcmp(semiology{i},'strong')
        semiology_coefficient(i) = 0.8489;
    else
        error('why')
    end

    % ictal coefficient
    if isempty(ictal{i})
        ictal_coefficient(i) = nan;
    elseif strcmp(ictal{i},'none/multilobar/diffuse')
        ictal_coefficient(i) = -0.8124;
    elseif strcmp(ictal{i},'focal')
        ictal_coefficient(i) = 0.8442;
    end

   

end

outT = table(names,lesion,interictal,neuropsych,semiology,ictal,final_focal,...
    intercept,lesion_coefficient,interictal_coefficient,neuropsych_coefficient,...
    semiology_coefficient,ictal_coefficient);

%% Do imputation
% restrict outT to allowable names
if impute_missing_data
    allowable_names = which_pts;
    allowable_rows = ismember(outT.names,allowable_names);
    outT(~allowable_rows,:) = [];
    
    outT.neuropsych_coefficient(isnan(outT.neuropsych_coefficient)) = ...
        nanmean(outT.neuropsych_coefficient);
end



 % Add the coefficients plus the intercept
sum_coefficients = outT.intercept + outT.lesion_coefficient + outT.interictal_coefficient + ...
    outT.neuropsych_coefficient + outT.semiology_coefficient + outT.ictal_coefficient;

% calculate the probability
prob = exp(sum_coefficients)./(exp(sum_coefficients)+1)*100;
outT.prob = prob;


% I confirmed that the 5 sense score is higher for the focal patients, nice
% positive control
if 0
    no_focal = cellfun(@isempty,final_focal);
    figure
    unpaired_plot(prob(strcmp(final_focal,'unifocal')),prob(strcmp(final_focal,'other')),...
        {'Focal','Non-focal'},'5-Sense score','para')
end
    
end


