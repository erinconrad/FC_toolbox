function outcome_num = parse_outcome_num(outcome,type)

if isempty(outcome)
    outcome_num = nan;
    return
end

if isnumeric(outcome) && isnan(outcome)
    outcome_num = nan;
    return
end

switch type
    case 'ilae'
        if contains(outcome,'1')
            outcome_num = 1;
        elseif contains(outcome,'2')
            outcome_num = 2;
        elseif contains(outcome,'3')
            outcome_num = 3;
        elseif contains(outcome,'4')
            outcome_num = 4;
        elseif contains(outcome,'5')
            outcome_num = 5;
        else
            outcome_num = nan;
        end
    case 'engel'
        
        if contains(outcome,'V')
            outcome_num = 4;
        elseif contains(outcome,'III')
            outcome_num = 3;
        elseif contains(outcome,'II')
            outcome_num = 2;
        elseif contains(outcome,'I')
            outcome_num = 1;
        else
            outcome_num = nan;
        end

end


end