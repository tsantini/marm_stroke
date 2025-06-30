%% demo for determination of a compression matrix from a calibration scan, followed by computation of trajectory from raw field probe data using compressed basis functions
%
%  This demo provides an example of computation of spherical harmonic
%  representation of phase over time based on measurements from field
%  probes. 
%
%  A standard least squares solution that also considers concomitant
%  gradient fields can be performed (see Vannesjo et al,
%  https://doi.org/10.1002/mrm.25841), or an algorithm that incorporates compressed basis functions in the fitting approach (manuscript in
%  preparation)
%
%  (c) Paul Dubovan, Corey Baron, 2024
%

% Load field probe data collected on the 7 Tesla MRI at the CFMM for a
% diffusion weighted calibration scan.
% Description of variables:
%   probePositions: position of each of the probes [m]
%   dataRaw:        time-series raw data measured by each field probe. The
%                   spherical harmonics are computed from its phase. 
%   dt:             time spacing between each raw data point [s]
%   fieldOffsets:   local field offset at each probe [T]. Measured using a
%                   calibration scan before the scan of interest
%   gammaProbes:    gyromagnetic ratio for probes. These ones used fluorine
%   gammaMRI:       gyromagnetic ratio for sample (typically proton)
%   B0:             magnetic field [T]
load data_calibrationscan_spiral.mat
probePositions = probePositions_100probes;
dataRaw = dataRaw_100probes_spiralpgse1p3mm_axial;
fieldOffsets = fieldOffsets_100probes;
gammaProbes = gammaProbes_100probes;
gammaMRI = gammaMRI_100probes;

%%
% Only keep weighted least squares correction on from previously proposed
% correction algorithm (Dubovan et al. https://doi.org/10.1002/mrm.29781)
opt.fitDoSteps = false;    % Stepwise correction
opt.fitDoWeights = true;  % Weighted least square correction (based on probe distance from isocenter)
opt.fitDoResAdj = false;   % Weighted residual correction

%% Ground truth fit: 100 probes
% Perform a fifth order fit using 100-probe array compiled from a
% calibration scan
opt.fitOrder = 5;  % Fit to 5th order spherical harmonics
[phs_spha_GT,phs_conc_GT] = harmonicsFromRaw(probePositions,dataRaw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt);

%% Compression Matrix Determination
% Define how many terms to include in the principal component analysis.
% For use of all terms calculated, keep set to 5
maxOrder = 5;
if maxOrder == 5
    mI = 36;
elseif maxOrder == 4
    mI = 25;
elseif maxOrder == 3
    mI = 16;
end
phs_spha = phs_spha_GT';
phs_spha = phs_spha(1:mI,:);

% Leave first order terms out for compression
phs_spha_new = phs_spha(5:mI,:);

% Scale each order by respective scaling factor
phs_spha_new(1:5,:) = phs_spha_new(1:5,:)/100; 
phs_spha_new(6:12,:) = phs_spha_new(6:12,:)/1000;
if maxOrder > 3
    phs_spha_new(13:21,:) = phs_spha_new(13:21,:)/10000;
end
if maxOrder > 4
    phs_spha_new(22:32,:) = phs_spha_new(22:32,:)/100000;
end

% Perform SVD
[compMat_temp, sv] = findCoilCompressMat(phs_spha_new,0,0);    

% Define a singular value threshold, and keep only the singular values up
% to the threshold
Nkeep = 5;
compMat = compMat_temp(1:Nkeep,:);

%%
% Load concurrent field probe data collected on the 7 Tesla MRI at the CFMM for a
% diffusion weighted in vivo scan.
load data_invivoscan_spiral.mat
probePositions = probePositions_16probes;
dataRaw = dataRaw_16probes_spiralpgse1p3mm_axial;
fieldOffsets = fieldOffsets_16probes;
gammaProbes = gammaProbes_16probes;
gammaMRI = gammaMRI_16probes;

%% Conventional fit: 16 probes
opt.fitOrder = 2;  % Fit to 2nd order spherical harmonics
% Perform second order fitting using conventional least sqaures method
[phs_spha_conventional, phs_conc_conventional] = harmonicsFromRaw(probePositions,dataRaw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt);

%% Compressed fit: 16 probes
opt.compMat = compMat; % Provide compression matrix to fitting algorithm
% Perform fit using compressed fitting method
[phs_spha_compressed, phs_conc_compressed] = harmonicsFromRaw(probePositions,dataRaw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt);

%%
% Plotting
figure;
NP = size(phs_spha_GT,2);
nr = floor(sqrt(NP));
nc = ceil(NP/nr);
t = dt*(1:size(phs_spha_GT,1))*1000; % [ms]
ylabelsAll = {'rad', 'rad/m', 'rad/m^2', 'rad/m^3', 'rad/m^4', 'rad/m^5'};

% Initialize empty arrays for plot handles
hGT = [];
hCompressed = [];
hConventional = [];

% Loop through and plot all k-coefficient terms
for np = 1:NP
    subplot(nr,nc,np);
    
    % Plot ground truth and store the handle
    hGT = plot(t,phs_spha_GT(:,np)); 
    hold('all');
    
    % Plot compressed data and store the handle
    hCompressed = plot(t,phs_spha_compressed(:,np));
    
    % Plot conventional data only for np < 10, and store the handle
    if np < 10
        hConventional = plot(t,phs_spha_conventional(:,np));
    end
    
    xlabel('time (ms)');
    
    % Set y-labels based on the value of np
    if np > 25
        ylabel(ylabelsAll{6});
    elseif np > 16
        ylabel(ylabelsAll{5});
    elseif np > 9
        ylabel(ylabelsAll{4});
    elseif np > 4
        ylabel(ylabelsAll{3});
    elseif np > 1
        ylabel(ylabelsAll{2});
    else 
        ylabel(ylabelsAll{1});
    end
end

legend([hGT, hCompressed, hConventional], {'Ground truth fit: 100 probes', 'Compressed fit: 16 probes', 'Conventional fit: 16 probes'});

  




