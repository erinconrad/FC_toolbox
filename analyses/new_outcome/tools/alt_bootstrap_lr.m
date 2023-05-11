function [conc_real,P,h] = alt_bootstrap_lr(T,features,nb)

% This is more bootstrapping for confidence intervals. I am trying to
% generate confidence intervals for the top feature sets for L and for R,
% and see how often the overlap is less than a chance overlap. 

rng(0) % seed a random number generator

% Establish output variables
nfeatures = length(features);
vec1 = cell(nfeatures,nb);
vec2 = cell(nfeatures,nb);
conc = nan(nfeatures,nb);

response = T.soz_lats;
npts = size(T,1);

%% Get the bootstrap best features
% loop over nb
for ib = 1:nb

    if mod(ib,100) == 0
        fprintf('\nDoing %d of %d\n',ib,nb);
    end

    % sample from the data with replacement (keeping original L/R/B
    % designations)
    which_pts = randi(npts,npts,1); % npt x 1 array drawn from 1:npts
    fakeT = T(which_pts,:);

    all_effects = nan(nfeatures,2);

    % loop over features
    for i = 1:nfeatures
        curr_feat = fakeT.(features{i});

        % Get effect sizes
        dL = cohenD(curr_feat(strcmp(response,'left')),curr_feat(strcmp(response,'bilateral')));
        dR = cohenD(curr_feat(strcmp(response,'right')),curr_feat(strcmp(response,'bilateral')));
        all_effects(i,:) = [abs(dL) abs(dR)]; % absolute value of effect size

    end

    % Ensure no nans
    assert(sum(isnan(all_effects),'all')==0)

    % Sort absolute value effect sizes in descending order
    [~,I1] = sort(all_effects(:,1),'descend');
    [~,I2] = sort(all_effects(:,2),'descend');

    % fill up feature names
    vec1(:,ib) = features(I1);
    vec2(:,ib) = features(I2);

    % get concordance at the top
    conc(:,ib) = concordance_at_the_top(vec1(:,ib),vec2(:,ib));

end

%% Do it once for real
real_effects = nan(nfeatures,2);
% loop over features
for i = 1:nfeatures
    curr_feat = T.(features{i});

    % Get effect sizes
    dL = cohenD(curr_feat(strcmp(response,'left')),curr_feat(strcmp(response,'bilateral')));
    dR = cohenD(curr_feat(strcmp(response,'right')),curr_feat(strcmp(response,'bilateral')));
    real_effects(i,:) = [abs(dL) abs(dR)]; % absolute value of effect size

end
% Sort absolute value effect sizes in descending order
[~,I1] = sort(all_effects(:,1),'descend');
[~,I2] = sort(all_effects(:,2),'descend');

% fill up feature names
vec1_real = features(I1);
vec2_real = features(I2);

% Get a CAT 
conc_real = concordance_at_the_top(vec1_real,vec2_real);

%% Get 95% CI for bootstrap conc
P = prctile(conc,[2.5 97.5],2);

%% Get a p-value - how often does 95% CI not include chance, in either direction?
chance_conc = ((1:nfeatures)/nfeatures)';
h = chance_conc < P(:,1) | chance_conc > P(:,2);

% correction for machine precision issues: if P upper and lower bounds are
% nearly the same, make h 0
h(abs(P(:,2)-P(:,1))<1e-3) = 0;

if 0
    figure
    shaded_error_bars_fc(1:nfeatures,conc_real,P,[])
    hold on
    plot(1:nfeatures,chance_conc,'k--','linewidth',2)

end


end