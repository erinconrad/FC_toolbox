library(pROC)

file_path = "/Users/erinconrad/Desktop/research/FC_toolbox/results/analysis/new_outcome/plots/delong_test.csv"
data1 <- read.csv(file_path)

# Get ROC
roc1 <- roc(data1$target,data1$pred_1)
roc2 <- roc(data1$target,data1$pred_2)

# do DeLong test
roc_comp <- roc.test(roc1,roc2,method="delong")

# when I do this for the ILAE outcomes, the AUCs for the combined and simple model are 0.857 and 0.691, respectively.
# DeLong test gives z = 1.7 and p = 0.088. This is the same as wha I get using the wilcoxonConfidence in Matlab.