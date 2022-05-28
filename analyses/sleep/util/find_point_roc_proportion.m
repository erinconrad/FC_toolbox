function closest_point = find_point_roc_proportion(X,Y,proportion)

%{
Proportion of positives = P/things
X = FPR = FP/N -> X*N = FP
Y = TPR = TP/P -> Y*P = TP
FP + TP = P
P + N = things


X*N + Y*P = P
X*(things-P) + Y*P = P
X*(1-P/things) + Y*P/things = P/things
X*(1-proportion) + Y*proportion = proportion


X*things - P*things + Y*P = P
X*things = P + P*things - Y*P
X*things = P(1+things-Y)
P = X*things/(1+things-Y)


X-X*proportion + Y*proportion = proportion
X = proportion + X*proportion - y*proportion
X = proportion(1+X-Y)
proportion = X/(1+X-Y)
%}

%% Find point on ROC curve that achieves desired proportion of positives
diff_to_min = abs(X*(1-proportion) + Y*proportion - proportion);
[~,closest_point] = min(diff_to_min);

%% Sanity test
test_prop = X(closest_point)/(1+X(closest_point)-Y(closest_point));
assert(abs(test_prop - proportion) < 0.005)

%% Get PPV and NPV
%{
NPV = TN/(TN + FN)
PPV = TP/(TP + FP)
PPV = (sensitivity x prevalence) / [ (sensitivity x prevalence) + ((1 – specificity) x (1 – prevalence)) ]
sensitivity = TPR
specificity = 1-FPR

NPV = (specificity x (1 – prevalence)) / [ (specificity x (1 – prevalence)) + ((1 – sensitivity) x prevalence) ]

PPV = TPR * proportion/[TPR * proportion + FPR *(1-proportion)]
NPV = (1-FPR)*(1-proportion)/((1-FPR)*(1-proportion) + (1-TPR)*proportion)
%}
PPV = Y(closest_point)*proportion/(Y(closest_point)*proportion + X(closest_point)*(1-proportion));
NPV = (1-X(closest_point))*(1-proportion)/...
    ((1-X(closest_point))*(1-proportion) + (1-Y(closest_point))*proportion);

end