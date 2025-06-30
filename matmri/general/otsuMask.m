function thresh = otsuMask(im)
% Find threshold to mask background using otsu's method

% Find histogram
numBins = max(100, round(numel(im)/10));
numBins = min(numBins,numel(im));
[counts,edges] = histcounts(abs(im(:)),numBins);

% Find class probabilities
w0 = cumsum(counts);
w1 = sum(counts) - w0;

% Find class means
tmp1 = (1:length(w0)).*counts;
tmp2 = cumsum(tmp1);
u0 = tmp2./w0;
u1 = (sum(tmp1) - tmp2)./w1;

% Find inter-class variance
sig2 = w0.*w1.*(u0-u1).^2;

% Find the threshold
[~,ind] = max(sig2);
thresh = edges(ind);

end