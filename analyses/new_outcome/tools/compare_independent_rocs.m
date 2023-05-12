function [p,z] = compare_independent_rocs(auc1,auc2,class1,class2,pos1,pos2)

% Following Mcclish, Med Decis Making 7:149-155, 1987, also, more famously
% Hanley and McNeil, The meaning and use of the area under a Receiver 
% Operating Characteristic (ROC) curve. Radiology, 1982, 143, 29-36.


npos1 = sum(strcmp(class1,pos1));
npos2 = sum(strcmp(class2,pos2));

nneg1 = sum(~strcmp(class1,pos1));
nneg2 = sum(~strcmp(class2,pos2));

calcQ1 = @(x) x/(2-x);
calcQ2 = @(x) 2*x^2/(1+x);
calc_var = @(x,n1,n2) (x*(1-x) + (n1-1)*(calcQ1(x)-x^2) + (n2-1)*(calcQ2(x)-x^2))/(n1*n2);

z = (auc1-auc2)/sqrt(calc_var(auc1,npos1,nneg1) + calc_var(auc2,npos2,nneg2));
p = 2*(1-normcdf(abs(z)));

end