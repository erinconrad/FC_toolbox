function prob = prob_fcn_for_carlos(L,R)

% establish coef (this will change). There will actually be two different
% models, each with a different coef
coef = [-0.4564 2.9101];

% Calculate AI
AI = (L-R)/(L+R);

% Do logistic function
prob = 1/(1+exp(-coef(1)+coef(2)*AI));

end