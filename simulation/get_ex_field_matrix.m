function [XImage,YImage,FieldMatrix] = ...
  get_ex_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters)
% This function calculates the field matrix A_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component, as well as the derivatives of A_{jk} w.r.t. the xyz coordinates
% of the emitter and w.r.t. the emission wavelength lambda.
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%
% parameters: NA, refractive indices of medium, wavelength (in nm), 
% nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction

% 2018 07 09, bug fix, MSie.

NA = parameters.NA;
lambda = parameters.lambda;
xemit = parameters.xemit;
yemit = parameters.yemit;
zemit = parameters.zemit;
xrange = parameters.xrange;
yrange = parameters.yrange;
zmin = parameters.zrange(1);
zmax = parameters.zrange(2);
ImageSizez = (zmax-zmin)/2;
Npupil = parameters.Npupil;
Mx = parameters.Mx;
My = parameters.My;
Mz = parameters.Mz;

% pupil and image size (in diffraction units)
PupilSize = 1.0;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

% image coordinate sampling (in physical length units).
DxImage = 2*ImageSizex/Mx;
DyImage = 2*ImageSizey/My;
ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
[YImage,XImage] = meshgrid(yimagelin,ximagelin);
XImage = (lambda/NA)*XImage;
YImage = (lambda/NA)*YImage;

% enlarging ROI in order to accommodate convolution with bead in
% computation PSF and PSFderivatives
XImage_tmp = XImage;
YImage_tmp = YImage;

% calculate auxiliary vectors for chirpz
[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,Mx);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,My);

% calculation Zernike mode normalization
if strcmp(parameters.fitmodel,'aberrations') 
  orders = parameters.aberrations(:,1:2);
  normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
  %normfac=1;
end

if strcmp(parameters.fitmodel,'aberrationsamp') 
  orders = parameters.aberrations(:,1:2);
  normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
  %normfac=1.0;
end

% set number of relevant parameter derivatives
switch parameters.fitmodel
  case 'xy'
    numders = 2;
  case 'xyz'
    numders = 3;
  case 'xylambda'
    numders = 3;
  case 'xyzlambda'
    numders = 4;
  case 'aberrations'
    numzers = size(parameters.aberrations,1);
    numders = 3+numzers;
  case 'aberrationsamp'
    numzers = size(parameters.aberrations,1);
    numders = 3+numzers-3+numzers;
end

% loop over emitter z-position
if Mz==1
  ZImage = (zmin+zmax)/2;
else
  DzImage = 2*ImageSizez/(Mz-1);
% enlarging axial range in order to accommodate convolution with bead
ZImage = linspace(zmin,zmax,Mz);
FieldMatrix = cell(3,2,Mz);


for jz = 1:numel(ZImage)
  zemitrun = ZImage(jz); 

% phase contribution due to position of the emitter
% Bug fix: minus of zemitrun replaced by plus sign and compex conjugate is 
% needed for wavevector{3}. MSie, 09 07 2018.

%   wavevector{3} = conj(wavevector{3});
  Wpos = xemit*wavevector{1}+yemit*wavevector{2}+zemit*wavevector{3};
  if strcmp(parameters.ztype,'medium')
    Wpos = Wpos+zemitrun*wavevector{3};
  end
  if strcmp(parameters.ztype,'stage')
    Wpos = Wpos+zemitrun*wavevectorzimm;
  end
  PositionPhaseMask = exp(-1i*Wpos);
%   PositionPhaseMask = exp(1i*Wpos);
  
  for itel = 1:3
    for jtel = 1:2
      % pupil functions and FT to matrix elements
      PupilFunction = PositionPhaseMask.*PupilMatrix{itel,jtel};
      IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
      FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
 
    end
  end
end

end