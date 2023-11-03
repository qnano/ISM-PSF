function [PSF] = get_ex_psfs_scalar(FieldMatrix,parameters)
% This function calculates the free or fixed dipole PSFs given the field
% matrix, the dipole orientation, and the pupil polarization, as well as
% the derivatives w.r.t. the xyz coordinates of the emitter and w.r.t. the
% emission wavelength lambda.
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%
% parameters: emitter/absorber dipole orientation (characterized by angles
% pola and azim), detection/illumination polarization in objective lens
% back aperture (characterized by angles alpha and beta).
pola = parameters.pola;
azim = parameters.azim;
polarizationpupil = parameters.polarizationpupil;
alpha = parameters.alpha;
beta = parameters.beta;

dipor(1) = sin(pola)*cos(azim);
dipor(2) = sin(pola)*sin(azim);
dipor(3) = cos(pola);

polpupil(1) = cos(alpha)*exp(1i*beta);
polpupil(2) = sin(alpha)*exp(-1i*beta);

% find dimensions and number of derivatives from input
dims = size(FieldMatrix);
if (length(dims)>3)
  Mz = dims(3);
  numders = dims(4);
  imdims = size(FieldMatrix{1,1,1});
elseif(length(dims)==3)
    Mz = dims(3);
    imdims = size(FieldMatrix{1,1,1});
else
  Mz = 1;
  numders = dims(3);
  imdims = size(FieldMatrix{1,1});
end
Mx = imdims(1);
My = imdims(2);
%disp(numders);
% calculation of free and fixed dipole PSF and the derivatives for the focal stack
PSF = zeros(Mx,My,Mz);


for jz = 1:Mz
% calculation of fixed PSF and derivatives

    Ex = FieldMatrix{1,1,jz};
    %Ey = FieldMatrix{2,2,jz};
    %Ez = FieldMatrix{3,1,jz}+FieldMatrix{3,2,jz};
  
    PSF(:,:,jz) = abs(Ex).^2;
    
    
  
end



% 3D convolution of the PSFs and derivatives with a bead


end

