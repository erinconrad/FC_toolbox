function [labels1,thing1,labels2,thing2] = find_matching_labels(labels1,thing1,labels2,thing2)

% find labels1 that are in labels2
[li12,loc12] = ismember(labels1,labels2);
labels2 = labels2(loc12(li12));
thing2 = thing2(loc12(li12));

% find labels2 that are in labels1
[li21,loc21] = ismember(labels2,labels1);
labels1 = labels1(loc21(li21));
thing1 = thing1(loc21(li21));

assert(isequal(labels1,labels2))


end