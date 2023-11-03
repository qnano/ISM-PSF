clear all
load('C:/code/vectorfitter/defult.mat');
sph=0;
n_pixels =31;
hn=(n_pixels-1)/2;
n_z_slices =3;
aberrations = [1,-1,0,0; 1,1,0.0,0; 2,0,0,0.0; 2,-2,0,0.0; 2,2,0.0,0.0];
count=1;
ma = 5;
for i=1:ma
    mind=-i:2:i;
    for j=1:length(mind)
        aberrations(count,:)=[i,mind(j),0,0];
        count=count+1;
    end
end
aberrations(3,2)=0;
aberrations(4,2)=-2;

aberrations(5,3)=72; % 72 nm astigmatism; aberration index follow the noll index third column in aberration stands for the phase aberration, and fourth column stands for the amplitude aberration


parameters.aberrations=aberrations;

parameters.NA = 1.49;



parameters.pixelsize=108.33; % nm
%parameters.samplingdistance = parameters.pixelsize;
parameters.fitmodel='aberrations';
%aberrations = [1,-1,0,0; 1,1,0.0,0; 2,0,0,0.05; 2,-2,0,0.0; 2,2,50.0,0.01; 3,-1,0.0,0.0; 3,1,0.0,0; 4,0,0.0,0.1; 3,-3,-0.0,0; 3,3,0.0,0; 4,-2,0.0,0; 4,2,0.0,0; 5,-1,0.0,0; 5,1,0.0,0; 6,0,0.0,0; 4,-4,0.0,0; 4,4,0.0,0;  5,-3,0.0,0; 5,3,0.0,0;  6,-2,0.0,0; 6,2,0.0,0; 7,1,0.0,0; 7,-1,0.0,0; 8,0,0.0,0];
parameters.beaddiameter = 40;
prameters.NA=1.33;
prameters.refimm=1.33;
parameters.lambda = 690;
parameters.lambdacentral=690;
parameters.lambdaspread=[690 690];
parameters.zrange = [-1500, 1500];
parameters.doetype='none';
parameters.dipoletype='free';
parameters.Mx = n_pixels ;
parameters.My = n_pixels ;
parameters.Mz = n_z_slices;
parameters.Npupil=41;
parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;

if strcmp(parameters.fitmodel,'aberrationsamp') 
    parameters.aberrations=aberrations;
    numparams_phaseaberr=length(aberrations)-3;
    numparams_ampaberr=length(aberrations);
    parameters.numparams=5+numparams_ampaberr+numparams_phaseaberr;
    parameters.numparams_phaseaberr=numparams_phaseaberr;
    parameters.numparams_ampaberr=numparams_ampaberr;
elseif strcmp(parameters.fitmodel,'aberrations') 
    aberrations=aberrations(4:end,1:3);
    parameters.aberrations=aberrations;
    parameters.numparams=5+length(aberrations);
end

for i=1:31
    for j=1:31
        parameters.xemit = (rand()*2-1)*108.33;
        parameters.yemit = (rand()*2-1)*108.33;
        [tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);


        [~,~,FieldMatrix,FieldMatrixDerivatives] = ...
          get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

        [PSF,PSFderivatives] = get_psfs_derivatives(FieldMatrix,FieldMatrixDerivatives,parameters);

    end
end




if false
    
    imgrow=(PSF.*1e6);
    img=zeros(n_pixels*2,n_pixels*2,n_z_slices);
    img(1:n_pixels,1:n_pixels,:)=imgrow;
    img(1:n_pixels,n_pixels+1:2*n_pixels,:)=imgrow;
    img(n_pixels+1:2*n_pixels,1:n_pixels,:)=imgrow;
    img(n_pixels+1:2*n_pixels,n_pixels+1:2*n_pixels,:)=imgrow;
    img=uint16(imgrow);
    zrange=parameters.zrange;
    fn="C:/issue1/supplymentary/sph_aberr_detection_rate/simulated_PSF.tif";
    %save(fn,'img','aberrations','zrange');
    for i=1:n_z_slices
        imwrite(img(:,:,i),fn,'WriteMode','append')
    end

end

