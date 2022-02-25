function plot_confusion_matrix(mat,class_names,xlab,ylab)

nclasses = length(class_names);

turn_nans_gray([1 0;0 1]);
colormap(gca,[0.8500, 0.3250, 0.0980;0.4660, 0.6740, 0.1880])
xticks(1:nclasses)
xticklabels(class_names)
yticks(1:nclasses)
yticklabels(class_names)
xlabel(xlab)
ylabel(ylab)
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d patients',mat(j,i)),'horizontalalignment','center','fontsize',20)
    end
end
set(gca,'fontsize',15);

end