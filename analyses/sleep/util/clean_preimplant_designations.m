function clean = clean_preimplant_designations(designation)

if contains(designation,'y','ignorecase',true)
    clean = 1;
elseif contains(designation,'n','ignorecase',true)
    clean = 0;
end

end