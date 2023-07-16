function multiclass_model_plots

%% Parameters
pca_all_perc = 95;
pca_spikes_perc = 95;
rm_non_temporal = 1;
rm_wake = 1;

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


%% Initialize figure
figure
set(gcf,'position',[1 1 1440 1000])
tiledlayout(2,3,"TileSpacing",'tight','padding','tight')

%% Run the lr_mt to extract features
if rm_wake == 1
    [T,features] =  lr_mt(3);
elseif rm_wake == 2
    [T,features] =  lr_mt(1);
else
    [T,features] =  lr_mt;
end
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% Identify HUP and MUSC
hup = contains(T.names,'HUP');
musc = contains(T.names,'MP');

%% Multiclass model - train on HUP and test on MUSC
just_spikes = 1;% all patients (spikes for now for speed)
combine_br = 0;
%hup = 
%out = classifier_wrapper(T,features,pca_all_perc,combine_br,just_spikes,rm_non_temporal,[],just_hup);
out = validation_classifier_wrapper(T,hup,musc,features,pca_all_perc,combine_br,just_spikes,rm_non_temporal);

%confusion matrix
nexttile
C = out.C;
classes = out.unique_classes;
nclasses = length(classes);
new_numbers = map_numbers_onto_range(C,[1 0]);
Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
D = diag(new_numbers);
Dcolor = [repmat(D,1,2),ones(length(D),1)];
Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
imagesc(Ccolor)
accuracy = sum(diag(C))/sum(C(:));
recall = nan(nclasses,1);
for i = 1:nclasses
    tp = C(i,i);
    fn = sum(C(i,~ismember(1:nclasses,i))); 
    recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
end
balanced_accuracy = mean(recall);

xticks(1:nclasses)
xticklabels((classes))
yticks(1:nclasses)
yticklabels((classes))
xlabel('Predicted')
ylabel('True')
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
    end
end

title(sprintf('Accuracy: %1.1f%%\nBalanced accuracy: %1.1f%%',...
    accuracy*100,balanced_accuracy*100))
set(gca,'fontsize',20)

end