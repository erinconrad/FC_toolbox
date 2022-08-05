function glm = outcome_logistic_regression(trainingData,forn)

switch forn
    case 'full'
        predictorNames = {'left other cortex proportion elecs', 'left temporal proportion elecs', 'right other cortex proportion elecs', 'right temporal proportion elecs', 'left other cortex proportion spikes', 'left temporal proportion spikes', 'right other cortex proportion spikes', 'right temporal proportion spikes'};
    case 'null'
        predictorNames = {'left other cortex proportion elecs', 'left temporal proportion elecs', 'right other cortex proportion elecs', 'right temporal proportion elecs'};
end

glm = fitglm(trainingData,'predictorvars',predictorNames,'responsevar','outcome','distribution','binomial');



end