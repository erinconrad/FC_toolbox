function [AI,score,pred] = online_calc_formula(left,right,coefs,classNames)

% Get AI
AI = (left-right)/(left+right);

% Get score
score = 1/(1+exp(-(AI*coefs(2)+coefs(1))));

% Get pred
if score >= 0.5
    pred = classNames{2};
elseif score < 0.5
    pred = classNames{1};
else
    error('what');
end

end