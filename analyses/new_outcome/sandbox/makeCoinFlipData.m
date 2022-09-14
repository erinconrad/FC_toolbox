function T = makeCoinFlipData(n_observations,n_predictors,predictor_type,response_type)


if strcmp(predictor_type,'binary') && strcmp(response_type,'binary')
    
    predictors = randi([0 1],n_observations,n_predictors);
    response = randi([0 1],n_observations,1);

    all = [response,predictors];
    T = array2table(all);
    T = renamevars(T,'all1','response');

end



end