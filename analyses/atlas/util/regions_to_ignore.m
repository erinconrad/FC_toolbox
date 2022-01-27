function ignore = regions_to_ignore(atlas,names)

ignore = zeros(length(names),1);

switch atlas
    case 'aal_bernabei'
        for i = 1:length(names)
            curr = names{i};
            
            % ignore if cerebelum or vernis
            if contains(curr,'Cerebelum'), ignore(i) = 1; end
            if contains(curr,'Vermis'), ignore(i) = 1; end
            
            % ignore if white matter
            if contains(curr,'White_Matter'), ignore(i) = 1; end
            
        end
        
end

ignore = logical(ignore);

end