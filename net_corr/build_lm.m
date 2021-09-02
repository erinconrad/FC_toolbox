function build_lm(all)

%{
Idea: what predicts the response at electrode j when you stimulate
electrode i? Some combination of the distance between them, the functional
connectivity between them, if one is a SOZ electrode or not


Switch things excluded for being subthreshold to zero??
%}
do_bin = 1;
which_tri = 'full';

%% Get three networks
pc = all.net.pc.data;
ccep = all.net.ccep.data;
dist = all.net.dist.data;

%% Make ccep nans zero
if do_bin
    is_nan = isnan(ccep);
    ccep(is_nan) = 0;
    ccep(~is_nan) = 1;
else
    ccep(isnan(ccep)) = 0;
end

%% Get outdegree, indegree, ns
switch which_tri
    case 'ut'
        tri = logical(triu(ones(size(ccep))));
     
    case 'lt'  
        tri = logical(tril(ones(size(ccep))));
    case 'full'
        tri = (ones(size(ccep)));
        tri(logical(eye(size(ccep)))) = 0;
        tri = tri(:);
        tri = logical(tri);
end

ccep_tri = ccep(tri);
dist_tri = dist(tri);
pc_tri = pc(tri);
        
tbl = table(ccep_tri,dist_tri,pc_tri);

if do_bin
    % Binary logistic regression model - is there a ccep response or not?
    mdl = fitglm(tbl,'ccep_tri ~ dist_tri + pc_tri','Distribution','binomial')
else
    % Linear model - how strong is the ccep response?
    mdl = fitlm(tbl,'ccep_tri ~ dist_tri + pc_tri')
end

end