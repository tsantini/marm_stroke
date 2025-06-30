**Package for MRI. Supports:**
  - Reconstruction
    - non-Cartesian regridding
    - iterative SENSE (Cartesian or non-Cartesian)
    - iterative SENSE for higher order models that may include: 
      - a B0 map
      - time varying spherical harmonics of phase accrual
  - diffusion MRI fitting
    - spatially regularized diffusion kurtosis fitting with an axially symmetric model (nii2kurt.m)
    - free water corrected kurtosis and micro-FA (nii2uFA_fwe.m). Cite with doi: 10.3389/fnins.2023.1074730
    
**Before Usage**
  - run setPath.m to add all the directories to the Matlab path
    
**Includes a basic image viewer and ROI drawing tool called "bview"**
  - enter "bview" in the Matlab command prompt to use
    
**Tips**
  - step through the demos in the demos folder
    - recommended order is demo_regridding, demo_regSENSE, demo_highOrder
  - open files to read usage info

**Notes**
  - Matlab gpuArray functionality is used whenever possible, so a compatible GPU is strongly recommended

**Acknowledgement**
  - Cite as: 
    1. Varela-Mattatall G, Dubovan PI, Santini T, Gilbert KM, Menon RS, Baron CA. Single-shot spiral diffusion-weighted imaging at 7T using expanded encoding with compressed sensing. Magn Reson Med. 2023 Apr 10. doi: 10.1002/mrm.29666
    2. Baron CA (2021, February 2). MatMRI: A GPU enabled package for model based MRI image reconstruction. Zenodo. http://doi.org/10.5281/zenodo.4495476
  - Additionally, please reference the following works for usage of the below functions:
    - **nii2kurt**: Hamilton, J., Xu, K., Geremia, N., Prado, V. F., Prado, M. A. M., Brown, A., & Baron, C. A. (2024). Robust frequency-dependent diffusional kurtosis computation using an efficient direction scheme, axisymmetric modelling, and spatial regularization. Imaging Neuroscience, 2, 1–22.
    - **nii2uFA_fwe**: Arezza NJJ, Santini T, Omer M, Baron CA. Estimation of free water-corrected microscopic fractional anisotropy. Front Neurosci. 2023 Mar 7;17:1074730. doi: 10.3389/fnins.2023.1074730. 
    - **harmonicsFromRaw**: Dubovan PI, Gilbert KM, Baron CA. A correction algorithm for improved magnetic field monitoring with distal field probes. Magn Reson Med. 2023 Dec;90(6):2242-2260. doi: 10.1002/mrm.29781
    - **findDelAuto**: Dubovan PI, Baron CA. Model-based determination of the synchronization delay between MRI and trajectory data. Magn Reson Med. 2022 Sep 26. doi:10.1002/mrm.29460.  
    - **nufftOp**: Baron CA, Dwork N, Pauly JM, Nishimura DG. Rapid compressed sensing reconstruction of 3D non-Cartesian MRI. Magn. Reson. Med. 2018;79:2685–2692.
    - **sampHighOrder**: 
      - Wilm BJ, Barmet C, Pruessmann KP. Fast higher-order MR image reconstruction using singular-vector separation. IEEE Trans Med Imaging. 2012 Jul;31(7):1396-403. doi: 10.1109/TMI.2012.2190991. Epub 2012 Mar 14. Erratum in: IEEE Trans Med Imaging. 2012 Sep;31(9):1833.
      - Baron CA, Dwork N, Pauly JM, Nishimura DG. Rapid compressed sensing reconstruction of 3D non-Cartesian MRI. Magn. Reson. Med. 2018;79:2685–2692.
    

(c) 2020, Corey Baron
