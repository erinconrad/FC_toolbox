function [f,vec_dens,vec_conn] = dens_model_specific_sr(all_locs,all_fc,all_soz_bin,all_spikes,max_spikes,sr)

%% establish all density matrices
npts = length(all_locs);
all_density = cell(npts,1);
for ip = 1:npts
    locs = all_locs{ip};
 
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % Make density matrix
    dens = interdistance_to_density_matrix(D,sr);
    all_density{ip} = dens;
    
end

[f,vec_dens,vec_conn] = inner_dens_model(all_density,all_fc,all_soz_bin,all_spikes,max_spikes);



end