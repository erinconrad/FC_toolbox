%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mni2cor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mni = cor2mni(cor, T)

if isempty(cor)
    mni = [];
    return;
end


mni = [cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]*((T'));
%coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]';
mni(:,4) = [];
%mni = round(mni);
return;