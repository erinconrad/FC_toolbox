function [engel1,engel2,ilae1,ilae2] = replace_with_my_outcomes(names,engel1,ilae1,engel2,ilae2,T)

npts = length(names);
table_names = T.Name;
for i = 1:npts
    % find matching row in table
    r = find(strcmp(names{i},table_names));

    if isempty(r), continue; end

    % replace engel and ilae
    engel1(i) = T.x1_yearEngel(r);
    engel2(i) = T.x2_yearEngel(r);

    ilae1{i} = sprintf('ILAE %d',T.x1_yearILAE(r));
    ilae2{i} = sprintf('ILAE %d',T.x2_yearILAE(r));


end


end