function show_confusion(C,class_labels,xlab,ylab,tit,acc)

%% Convert to proportion of total in predicted class
C = C./sum(C,1);

%% Define colors
nclasses = length(class_labels);

% Map numbers onto 0 to 1
new_numbers1 = map_numbers_onto_range(C,[0 1]);
new_numbers2 = map_numbers_onto_range(C,[1 0]);

%% First assign everything to be red
Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers2,1,1,2));

%% Next, find diagonal elements and make them blue
D = diag(new_numbers2);
Dcolor = [repmat(D,1,2),ones(length(D),1)];
Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;

imagesc(Ccolor)
hold on
for ic = 1:nclasses
    for jc = 1:nclasses
        text(ic,jc,sprintf('%1.2f',C(jc,ic)),...
            'horizontalalignment','center','fontsize',20,'fontweight','bold')
    end
end
xticks(1:nclasses)
yticks(1:nclasses)
xticklabels(class_labels)
yticklabels(class_labels)
xlabel(xlab)
ylabel(ylab)
set(gca,'fontsize',20)
title(sprintf('%s Average accuracy: %1.1f%%',tit,mean(acc)*100))



end