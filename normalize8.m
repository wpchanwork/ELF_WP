function I8 = normalize8(I)
    % Convert to uint8 safely (handles HDR / out-of-[0,1] ranges)
    if isa(I,'uint8'), I8 = I; return; end
    I = double(I);
    if all(I(:)>=0 & I(:)<=1)
        I8 = im2uint8(I);
    else
        lo = prctile(I(:),1); hi = prctile(I(:),99);
        I  = (I - lo) / max(hi - lo, eps);
        I8 = im2uint8(min(max(I,0),1));
    end
end