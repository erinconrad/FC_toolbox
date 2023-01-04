function ll = line_length(x)

y = x(1:end-1,:);
z = x(2:end,:);

ll = sum(abs(z-y),1)/size(x,1);


end