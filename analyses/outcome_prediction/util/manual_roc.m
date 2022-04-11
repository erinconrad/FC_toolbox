function W = manual_roc(true,class)

% Only works for two classes.
% nN = number normal 
% nA = number abnormal  
N = find(true == 0);
A = find(true == 1);
nN = sum(true==0);
nA = sum(true==1);

% Build S matrix

S = nan(nA,nN);

% loop over all normals
for in = 1:nN
    
    % and all abnormals
    for ia = 1:nA
        
        xA = class(A(ia));
        xN = class(N(in));
        
        % compare 
        if xA > xN
            S(ia,in) = 1;
        elseif xA == xN
            S(ia,in) = 0.5;
        elseif xA < xN
            S(ia,in) = 0;
        end
        
    end

end

W = 1/(nN*nA) * sum(S,'all'); % AUC
%Q1 = W/(2-W);
%Q2 = 2*W^2/(1+W);

%SE = sqrt((W*(1-W)+(nA-1)*(Q-W^2)+(nN-1)*(Q2-W^2))/(nA*nN));


end