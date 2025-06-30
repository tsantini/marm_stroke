function [AtA,Aty] = mrSampFuncMat(y,Crcvr,b0,sampTimes,phs_spha,phs_conc,phs_grid,imMask)
% Create A'*A and A'*y, where A is the dense MRI sampling matrix and y is the data vector. 
%   Modern GPUs make it feasible to store the dense array for a single slice, 
%   making non-iterative solutions possible 
%
%  [AtA,Aty] = mrSampFuncMat(data,Crcvr,b0,sampTimes,phs_spha,phs_conc,phs_grid)
%
% Outputs:
%   AtA: dense matrix that represents C'*S'*S*C, where C applies reciever
%       sensitivities and S performs MRI sampling
%   Aty: vector output of C'*S'*data, where data is the k-space data
%
% Inputs:
%   data:  kspace data
%   Crcvr: receiver sensitivity profiles
%   b0,sampTimes,phs_spha,phs_conc,phs_grid: k-space and B0 inhomogeneity
%     data. See sampHighOrder documentation for more details.
%
% (c) Corey Baron 2022
%

% Create sampling operator that does not include receivers
% TODO: imMask could be an input here for a slight speedup
S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid);

% Build receiver operator
R = rcvrOp(Crcvr,0);

% Create data vector output
Aty = R'*(S'*y);
Aty = Aty(imMask(:));

% Create R'*S'*S*R
R = single(R.maps);
S = reshape(single(S.kbase), numel(b0), []).';
S = S(:,imMask(:));
S = exp(1i*S)/sqrt(numel(Aty));
S = S'*S;
AtA = 0;
for n=1:size(R,4)
    tmp = R(:,:,:,n);
    tmp = reshape(tmp(imMask(:)), 1, []);
    AtA = AtA + tmp'.*S.*tmp; 
end



end