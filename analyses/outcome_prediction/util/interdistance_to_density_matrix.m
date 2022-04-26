function dens = interdistance_to_density_matrix(dist,sr)
W = 1;
nchs = size(dist,1);
dens = nan(size(dist));
for i = 1:nchs
    for j = 1:i-1
         d = dist(i,j); 
         if d > sr
             tdens = 0;
         else
             tdens = (3/pi*W*(1-(d/sr)^2)^2)/sr^2;
         end
         
         dens(i,j) = tdens;
         dens(j,i) = tdens;
                
    end
end

end