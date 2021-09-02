function x = upper_unwrap(y)
    [dx, dy] = size(y);
    if dy == 1
        nchs = round(0.5 + sqrt(0.25+2*dx));
    else
        nchs = dy;
    end
    [X, Y] = ndgrid(1:nchs, 1:nchs);
    f = X < Y;
    if dy == 1
        x = zeros(nchs, nchs);
        x(f) = y;
        x = x + x';
    else
        x = y(f);
    end
end
