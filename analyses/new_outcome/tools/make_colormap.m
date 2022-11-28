function cvars_map = make_colormap(n,cvars)

cvars_map = round((cvars - min(cvars,[],'all'))/(max(cvars,[],'all')-min(cvars,[],'all')) * (n-1) + 1);


end