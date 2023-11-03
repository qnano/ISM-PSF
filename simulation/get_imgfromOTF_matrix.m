function [PSF] = ...
  get_imgfromOTF_matrix(OTF,parameters)


NA = parameters.NA;
lambda = parameters.lambda;
Npupil = parameters.Npupil;
Mx = parameters.Mx;
My = parameters.My;


xrange = parameters.pixelsize*parameters.Mx/2;
yrange = parameters.pixelsize*parameters.My/2;

PupilSize = 4.0;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,Mx);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,My);

IntermediateImage = transpose(cztfunc(OTF,Ay,By,Dy));
PSF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
PSF=abs(PSF);
end