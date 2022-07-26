function outcome_num = parse_outcome(outcome,type)

% 1 = good, 0 = bad
if isempty(outcome)
    outcome_num = nan;
    return
end

switch type
    case 'ilae'
        if contains(outcome,'1')
            outcome_num = 1;
        else
            outcome_num = 0;
        end
            
        
end
        

end