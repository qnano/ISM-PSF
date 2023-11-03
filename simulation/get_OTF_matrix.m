function [OTF] = ...
  get_OTF_matrix(PSF,parameters)


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

[Ax,Bx,Dx] = prechirpz(ImageSizex,PupilSize,Mx,Npupil);
[Ay,By,Dy] = prechirpz(ImageSizey,PupilSize,My,Npupil);

IntermediateImage = transpose(cztfunc(PSF,Ay,By,Dy));
OTF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

end