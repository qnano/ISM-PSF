function [OTF] = ...
  get_OTF_drom_img(img,parameters)


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

[Ax,Bx,Dx] = prechirpz(ImageSizex,PupilSize,Mx,Npupil);
[Ay,By,Dy] = prechirpz(ImageSizey,PupilSize,My,Npupil);

IntermediateImage = transpose(cztfunc(img,Ay,By,Dy));
OTF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

end