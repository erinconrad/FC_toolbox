function test_online_calc

%% Parameters
which_refs = {'machine','car','bipolar'};

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load mt_data
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% get variables of interest
names = mt_data.all_names;
spikes = mt_data.all_spikes(:,:,3);
labels = mt_data.all_labels;

% Loop over references
for ir = 1:length(which_refs)

    ref_spikes = spikes(:,ir);
    ref_labels = labels(:,ir);
    
    % Load model output file
    out = load([plot_folder,sprintf('ext_models_%s.mat',which_refs{ir})]);
    out = out.all;
    model = out.approach(1).model;
    model2 = model(2).val(2); 

    % Loop over sides
    for is = 1:2

        curr = model2.side(is).result;

        % Get the model info
        model_names = curr.names;
        model_scores = curr.scores;
        model_preds = curr.all_pred;
        model_coefs = curr.tc.coef;
        model_classNames = curr.unique_classes;
        model_classes = curr.class;

        % Initialize my test things
        test_scores = nan(length(model_scores),1);
        test_preds = cell(length(model_preds),1);
        test_left = nan(length(model_scores),1);
        test_right = nan(length(model_scores),1);
        test_ai = nan(length(model_scores),1);

        % Loop over the patients
        for ip = 1:length(model_names)
            
            % Find matching name in mt_data
            im = strcmp(model_names{ip},names);

            % get pt specific spikes and labels
            pt_spikes = ref_spikes{im};
            pt_labels = ref_labels{im};
            
            % Get left and right spike rates and AI
            %{
            [ai,left,right] = calc_ai_ns(pt_labels,pt_spikes,[],[],[],[],1,1,...
                'spikes',[],0);
            %}

            [ai,left,right] = clean_ai_calc(pt_labels,pt_spikes,1,1);

            % Get new score and pred
            [AI,score,pred] = online_calc_formula(left,right,model_coefs,model_classNames);

            % Make sure AIs agree
            assert(abs(AI-ai)<1e-4)

            test_scores(ip) = score;
            test_preds(ip) = {pred};
            test_left(ip) = left;
            test_right(ip) = right;
            test_ai(ip) = AI;

        end

        % Compare scores and preds
        assert(all(abs(model_scores-test_scores)<1e-3))
        assert(isequal(model_preds,test_preds))

        if 1
            fprintf('\n%s ref, side %d',which_refs{ir},is)
            table(model_names,model_classes,...
                cellfun(@(x,y) isequal(x,y),model_preds,model_classes),...
                test_left,test_right,test_ai,...
                test_scores,model_scores,test_preds,model_preds)
            curr.C
            curr.balanced_accuracy
        end
    
    end

end


end