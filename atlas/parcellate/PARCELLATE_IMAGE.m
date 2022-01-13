function X = PARCELLATE_IMAGE(target,template,nparc)

% INPUTS:
% target: X-by-Y-by-Z image
% template: X-by-Y-by-Z matrix of indices for parcels or groups within same
% space as target
% nparc: number of nodes in template to parcellate. will use ROIs labeled 1:nparc 

% ACTIONS:
% parcellate target 3D image according to indices in template

% OUTPUTS:
% returns vector X of mean values in each unique index of template image

if ~exist('nparc','var')
    nparc = max(unique(template));
end

X = zeros(nparc,1);
for ROI = 1:nparc
    X(ROI) = mean(mean(mean(target(template == ROI))));
end
