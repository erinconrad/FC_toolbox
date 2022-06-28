function out = parse_outcomes(outcome,type)

if isempty(outcome)
    out = nan;
    return
end

switch type
    case 'engel'
        
        
    case 'ilae'
        assert(strcmp(outcome(1:5),'ILAE '))
        
        
end

end