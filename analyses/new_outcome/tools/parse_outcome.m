function [outcome_num,rule] = parse_outcome(outcome,type)

% 1 = good, 0 = bad


switch type
    case 'ilae'
        if isempty(outcome)
            outcome_num = nan;
        elseif contains(outcome,'1') || contains(outcome,'2')
            outcome_num = 1;
        else
            outcome_num = 0;
        end
        rule = 'ILAE score of 1 or 2';
        
    case 'engel'
        if isempty(outcome)
            outcome_num = nan;
        elseif strcmp(outcome,'IA') || strcmp(outcome,'IB') || ...
                strcmp(outcome,'IC') || strcmp(outcome, 'ID')
        
            outcome_num = 1;
        else
            outcome_num = 0;
        end
        rule = 'Engel score IA-ID';  
        
end
        

end