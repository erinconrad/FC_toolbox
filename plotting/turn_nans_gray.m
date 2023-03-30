function h = turn_nans_gray(im,special_color,target)

    if exist('special_color','var') == 0
        special_color = [0.7 0.7 0.7];
    end


        if exist('target','var') == 0
            target = gca;
        end

    cmap = colormap(target);
    h = imagesc(im);
    % white
    
    nanjet = [special_color; cmap  ];
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
    
    set(target,'colormap',nanjet);
    if cmin == cmax && cmin == 1
        cmin = 0;
    end
    caxis([cmin cmax]);
end