function [rankOut, bvecFromMat, bval] = bmatRank(bmat,eigThresh)
% Approximates the rank of a list of b-matrices. Dim 1 should be b-matrix
% entries (s/mm2), dim 2 should be number of acquisitions. 
% %  6 matrix entries makes symmetric matrices be assumed. First three are
% diagonal (xx, yy, zz), next are cross terms (xy, xz, yz).
%
% b = 0 scans, where all eigenvalues will be below the threshold, will have rank of 1

% Set threshold used for eigenvalues (s/mm^2)
if nargin<2
    % eigThresh of 15 means that the smallest STE b-shell detectable is
    % b=45 s/mm2
    eigThresh = 15;
end

rankOut = zeros(1,size(bmat,2),'like',bmat);
bvecFromMat = zeros(3,size(bmat,2),'like',bmat);
bval = zeros(1,size(bmat,2),'like',bmat);
for n=1:size(bmat,2)
    if size(bmat,1)==6
        % Bxx Byy Bzz Bxy Bxz Byz
        bmat_a = [bmat(1,n), bmat(4,n), bmat(5,n);
                  bmat(4,n), bmat(2,n), bmat(6,n);
                  bmat(5,n), bmat(6,n), bmat(3,n)];
    elseif size(bmat,1)==9
        % Bxx Bxy Bxz Bxy Byy Byz Bxz Byz Bzz
        bmat_a = [bmat(1,n), bmat(2,n), bmat(3,n);
                  bmat(4,n), bmat(5,n), bmat(6,n);
                  bmat(7,n), bmat(8,n), bmat(9,n)];
    else
        error('unkown size of dim 1')
    end

    % Find eigenvalues
    [v,e] = eig(bmat_a);
    bvecFromMat(:,n) = v(:,3);
    bval(n) = sum(diag(e));

    % Find rank
    rankOut(n) = sum(diag(e)>eigThresh);
    
    % Account for b=0 scans
    if rankOut(n) < 1
        rankOut(n) = 1;
    end
end