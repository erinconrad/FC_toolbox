function [outcome_num,rule] = parse_outcome_new(outcome,type)

% 1 = good, 0 = bad


switch type
    case 'ilae'
        if isempty(outcome)
            outcome_num = '';
        elseif isnumeric(outcome) && isnan(outcome)
            outcome_num = '';
        elseif contains(outcome,'1') || contains(outcome,'2')
            outcome_num = 'good';
        else
            outcome_num = 'bad';
        end
        rule = 'ILAE score of 1 or 2';
        
    case 'engel'
        if isempty(outcome)
            outcome_num = '';
        elseif isnumeric(outcome) && isnan(outcome)
            outcome_num = '';
            
        elseif strcmp(outcome,'IA') || strcmp(outcome,'IB') || ...
                strcmp(outcome,'IC') || strcmp(outcome, 'ID')
 
            outcome_num = 'good';
        else
            outcome_num = 'bad';
        end
        rule = 'Engel score IA-ID';  
        
end
        

end