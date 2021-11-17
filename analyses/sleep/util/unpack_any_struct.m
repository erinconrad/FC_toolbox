function unpack_any_struct(S)


fields = fieldnames(S);
for f = 1:length(fields)
    fname = fields{f};
    assignin('caller',fname,S.(fname));
end

end