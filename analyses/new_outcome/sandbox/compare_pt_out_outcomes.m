function compare_pt_out_outcomes(pt,out)

pt_names = cell(length(pt),1);
out_names = cell(length(out.all_names),1);
pt_ilae_two = cell(length(pt),1);
out_ilae_two = cell(length(out.all_names),1);

for i = 1:length(pt)
    pt_names{i} = pt(i).name;
    pt_ilae_two{i} = pt(i).clinical.ilae{2};
    
end

for i = 1:length(out.all_names)
    out_names{i} = out.all_names{i};
    out_ilae_two{i} = out.all_two_year_ilae{i};
end

% Reconcile
[ia,ic] = ismember(out_names,pt_names);

table(out_names,out_ilae_two,pt_names(ic(ia)),pt_ilae_two(ic(ia)))

end