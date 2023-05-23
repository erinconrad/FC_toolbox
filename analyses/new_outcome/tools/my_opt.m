function [T_opt,I] = my_opt(X,Y,T)

to_opt = (1-X)+Y;
[max_opt,I] = max(to_opt);
T_opt = T(I);

end