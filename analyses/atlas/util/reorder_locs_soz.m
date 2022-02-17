function out = reorder_locs_soz(A,soz)

%{
Goal: reorder the SOZ matrix to be in the form of SOZ-ipsi non SOZ 1-ipsi
non SOZ2-contra to SOZ- contra to 1-contra to 2

Pseudocode:
loop over patients
find soz
find the reordering for that patient
do the reordering for that apteitn
%}

code = {'L1','R1','L2','R2','L3','R3'};

out = nan(size(A));
all_orders = nan(size(soz));
npts = size(A,3);

for ip = 1:npts
    % get soz
    curr_soz = soz(:,ip);
    
    % skip if empty
    if sum(curr_soz) == 0, continue; end
    
    % initialize new order
    new_order = nan(6,1);
        
    switch code{curr_soz}
        
        case 'L1'
            new_code = {'L1','L2','L3','R1','R2','R3'};
        case 'L2'
            new_code = {'L2','L3','L1','R2','R3','R1'};
        case 'L3'
            new_code = {'L3','L1','L2','R3','R1','R2'};
        case 'R1'
            new_code = {'R1','R2','R3','L1','L2','L3'};
        case 'R2'
            new_code = {'R2','R3','R1','L2','L3','L1'};
        case 'R3'
            new_code = {'R3','R1','R2','L3','L1','L2'};
            
    end
    
    for i = 1:length(new_order)
        
        new_order(i) = find(ismember(code,new_code{i}));
        
    end
    all_orders(:,ip) = new_order;
    
    out(:,:,ip) = A(new_order,new_order,ip);
    
end

if 0
    figure
    turn_nans_gray(nanmean(out,3))
    
end

end