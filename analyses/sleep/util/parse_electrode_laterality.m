function lat = parse_electrode_laterality(labels)

first_letter = cellfun(@(x) x(1),labels,'uniformoutput',false);
lat = cell(length(labels),1);
lat(strcmp(first_letter,'L')) = {'L'};
lat(strcmp(first_letter,'R')) = {'R'};


end