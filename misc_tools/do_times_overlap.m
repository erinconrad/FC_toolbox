function overlap = do_times_overlap(start_1,end_1,start_2,end_2)

% The only way they do not overlap is if end_1 is < start_2 (implying 1 is
% entirely before 2) or if start_1 is > end_2 (implying 1 is entirely after
% 2)

%% Double check that ends are after starts
assert(end_1>=start_1 && end_2>=start_2)

overlap = 1; % start assuming overlap

%% Is 1 before 2?
if end_1<start_2
    overlap = 0;
end

%% Is 1 after 2
if start_1 > end_2
    overlap = 0;
end

end