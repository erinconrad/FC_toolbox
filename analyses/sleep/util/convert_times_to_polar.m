function polar = convert_times_to_polar(times,type)

%{
I want midnight to be pi/2, 6 am to be 0, 12 pm to be -pi/2, 6 pm to be pi
%}
switch type
    case 'radians'
        polar = (times)/(24*3600)*2*pi;
        polar = polar;

    case 'degrees'
        
end
        

end