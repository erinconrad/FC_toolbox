function mdl = sozEnsemble(trainingData,forn)

inputTable = trainingData;
switch forn
    case 'full'
        predictorNames = {'left other cortex proportion elecs', 'left temporal proportion elecs', 'right other cortex proportion elecs', 'right temporal proportion elecs', 'left other cortex proportion spikes', 'left temporal proportion spikes', 'right other cortex proportion spikes', 'right temporal proportion spikes'};
    case 'null'
        predictorNames = {'left other cortex proportion elecs', 'left temporal proportion elecs', 'right other cortex proportion elecs', 'right temporal proportion elecs'};
end

predictors = inputTable(:, predictorNames);
response = inputTable.SOZ;

mdl = fitcensemble(...
    predictors, ...
    response, ...
    'ClassNames', {'bilateral/diffuse'; 'left other cortex'; 'left temporal'; 'right other cortex'; 'right temporal'});


end