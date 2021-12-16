function norm_thing = norm_exc(thing,thing_for_norm,exc)

exc_bin = zeros(length(thing_for_norm),1);
exc_bin(exc) = 1;
norm_thing = (thing - nanmedian(thing_for_norm(~exc_bin)))./iqr(thing_for_norm(~exc_bin));

end