function h = turn_nans_gray(im)
    % white
    cmap = colormap;
    nanjet = [ 0.7,0.7,0.7; cmap  ];
    nanjetLen = length(nanjet); 
    pctDataSlotStart = 2/nanjetLen;
    pctDataSlotEnd   = 1;
    pctCmRange = pctDataSlotEnd - pctDataSlotStart;

    dmin = nanmin(im(:));
    dmax = nanmax(im(:));
    dRange = dmax - dmin;   % data range, excluding NaN

    cLimRange = dRange / pctCmRange;
    cmin = dmin - (pctDataSlotStart * cLimRange);
    cmax = dmax;
    h = imagesc(im);
    set(gcf,'colormap',nanjet);
    caxis([cmin cmax]);
end