function show_native_mni_locs(summ)

labels = summ.labels;
mni_locs = summ.locs;
native_locs = summ.native_locs;

figure
nexttile
scatter3(mni_locs(:,1),mni_locs(:,2),mni_locs(:,3),'white');
hold on
text(mni_locs(:,1),mni_locs(:,2),mni_locs(:,3),labels);


nexttile
scatter3(native_locs(:,1),native_locs(:,2),native_locs(:,3),'white');
hold on
text(native_locs(:,1),native_locs(:,2),native_locs(:,3),labels);
end