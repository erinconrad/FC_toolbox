function new_numbers = map_numbers_onto_range(numbers,range)

maxN = max(numbers,[],'all');
%minN = min(numbers,[],'all');
minN = 0;

new_numbers = (numbers - minN)./(maxN - minN)*(range(2)-range(1))+range(1);

end