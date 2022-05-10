function pretty_exp = pretty_exp_html(num)

exp_part = floor(log(abs(num))/log(10));
num_part = num/(10^exp_part);

pretty_exp = sprintf('%1.1f x 10<sup>%d</sup>',num_part,exp_part);

end