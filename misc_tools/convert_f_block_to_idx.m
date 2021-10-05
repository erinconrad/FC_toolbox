function idx = convert_f_block_to_idx(pt,p,f,b)

% This is a tool that takes, for a given patient p, a file number, f, and a
% block b within that file, and it outputs idx, the overall index for the
% segment

idx = 0;

% Loop through preceding files
for ifile = 1:f-1
    idx = idx + size(pt(p).ieeg.file(ifile).block_times,1); % add the prior blocks
end

% now add the current file block
idx = idx + b;

end