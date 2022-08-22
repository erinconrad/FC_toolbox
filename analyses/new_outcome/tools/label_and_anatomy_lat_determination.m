function lat = label_and_anatomy_lat_determination(labels,anatomy)

nchs = length(labels);
lat = cell(nchs,1);
for ich = 1:nchs
    %% First try anatomy
    ana = anatomy{ich};
    if isempty(ana)
        ana = ' ';
    end
    if strcmpi(ana(1),'L')
            lat{ich} = 'L';
    elseif strcmpi(ana(1),'R')
        lat{ich} = 'R';
    else
        %% Try label
        lab = labels{ich};
        if strcmpi(lab(1),'L')
            lat{ich} = 'L';
        elseif strcmpi(lab(1),'R')
            lat{ich} = 'R';
        else
            lat{ich} = '';
        end
    end
    
    %{
    %% First, see if the label starts with R or L
    lab = labels{ich};
    if strcmpi(lab(1),'L')
        lat{ich} = 'L';
    elseif strcmpi(lab(1),'R')
        lat{ich} = 'R';
    else
        %% if label ambiguous, check anatomy
        ana = anatomy{ich};
        if isempty(ana)
            lat{ich} = '';
            continue
        end
        if strcmpi(ana(1),'L')
            lat{ich} = 'L';
        elseif strcmpi(ana(1),'R')
            lat{ich} = 'R';
        else
            lat{ich} = '';
        end
    end
    %}
        
    
end

end