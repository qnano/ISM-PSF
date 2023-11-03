function [PSF] = ...
  get_PSFfromOTF_matrix(OTF,parameters)


NA = parameters.NA;
lambda = parameters.FTlambda;
Npupil = parameters.Npupil;
Mx = parameters.Mx;
My = parameters.My;

xrange = parameters.xrange;
yrange = parameters.yrange;

PupilSize = 4.0;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,Mx);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,My);

IntermediateImage = transpose(cztfunc(OTF,Ay,By,Dy));
PSF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
PSF=abs(PSF);
end