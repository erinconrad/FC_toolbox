function [sleep,wake,norm_ad] = find_sleep_wake(ad,exc,disc)

norm_ad = norm_exc(ad,ad,exc);

norm_ad(exc) = nan;
sleep = norm_ad <= disc;
wake = norm_ad > disc;

end