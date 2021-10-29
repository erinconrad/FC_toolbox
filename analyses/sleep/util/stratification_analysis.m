function all_ps = stratification_analysis(thing)

% thing is ngroups x npts x 2 where 2 is soz vs not soz
all_ps = size(thing,1);

for s = 1:size(thing,1)
    curr = squeeze(thing(s,:,:));
    
    all_ps(s) = signrank(curr(:,1),curr(:,2));
    
end

end