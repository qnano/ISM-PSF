%% ExamplePhaseRetrieval
% This codes shows an example of a TFS of a bead (175nm, green) with the
% deformable mirror (MIRAO52E) upon start up. the phase retrieval is used
% to estimate the Zernike coefficients. Subsequently, the negative values
% of these coefficients are then modulated by the DM and a second TFS is
% acquired. Phase retrieval on the corrected TFS shows that the Zernike
% modes are corrected.
% At each z-position 3 frames are acquired, resulting in 3 TFS per time.

% Marijn Siemons, 28-02-2020

% Fits correspond well within z-positions of +-500nm. The corrected TFS
% shows some minor deviations at large z-stage positions, probably due to 
% Zernike modes which are present in the wavefront, but not fitted by the 
% model. Also, the z-stage position might by slightly different in the
% experimental TFS which can lead to model errors.
clear all
close all
maxn=5;
count=1;
for i=1:maxn
    mind=-i:2:i;
    for j=1:length(mind)
        aberrations(count,:)=[i,mind(j),0,0];
        count=count+1;
    end
end
aberrations(3,2)=0;
aberrations(4,2)=-2;
aberrationsall=aberrations;
aberrations_amp=aberrations;
aberrations_phase=aberrations;
%aberrations =matfile.Abcorrection;
if false
    load('C:/research/chisquareeval/data/Exp003_correctedaberrationsand60mlambdaadded_Greenbead_180nm.mat');
    parameters.aberrationcorrected=0;
    parameters.NA = 1.49;
    parameters.pixelsize=80;
    parameters.samplingdistance = parameters.pixelsize;
    parameters.beaddiameter = 180;
    parameters.lambda = 690;
    imall=zeros(31,31,NumberImages,5,'uint16');
    
    %% load
    %imgind=6;
    for imgind=6:6
        im=zeros(31,31,NumberImages,'uint16');
        if (imgind<10)
            %FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/0"+num2str(imgind)+".tif";
            FileTif="C:/code/ExamplePhaseRetrieval/data/simdata/testimg.tif";
            %FileTif="C:/research/chisquareeval/data/chivsI9exp_fl/zernike_mode_img/testimg10/roi301zernike2250_13.tif";
        else
            FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/"+num2str(imgind)+".tif";
        end
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        NumberImages=length(InfoImage);


        for i=1:NumberImages
           im(:,:,i)=imread(FileTif,'Index',i);
        end
        im=double(im);
        im=im./1e6;
        im=im.*5e3;
        im=poissrnd(im+3);
    end
    TFSuncorexp=im;
else
    
    load('C:/research/chisquareeval/data/Exp003_correctedaberrationsand60mlambdaadded_Greenbead_180nm.mat');
    parameters.aberrationcorrected=0;
    parameters.samplingdistance = parameters.pixelsize;

    
    imall=zeros(31,31,21,3,'uint16');
    %% load
    %imgind=6;
    for imgind=6:8
        if (imgind<10)
            FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/0"+num2str(imgind)+".tif";
            %FileTif="C:/code/ExamplePhaseRetrieval/data/simdata/testimg.tif";
            %FileTif="C:/research/chisquareeval/data/chivsI9exp_fl/zernike_mode_img/testimg10/roi301zernike2250_13.tif";
        else
            FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/"+num2str(imgind)+".tif";
        end
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        NumberImages=length(InfoImage);

        im=zeros(nImage,mImage,NumberImages,'uint16');
        for i=1:NumberImages
           im(:,:,i)=imread(FileTif,'Index',i);
        end
        im=double(im);

        im = (double(im) - cameraoffset).*cameragain;
        ROIsize = parameters.Mx;
        hroi=(ROIsize-1)/2;
        loc=[66 66];
        im=im(loc(1)-hroi:loc(1)+hroi,loc(2)-hroi:loc(2)+hroi,:);
        imall(:,:,:,imgind-5)=im;
    end
    
    TFSuncorexp=double(imall);

end

im=zeros(31,31,NumberImages,'uint16');
imgind=6;
if (imgind<10)
    %FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/0"+num2str(imgind)+".tif";
    FileTif="C:/code/ExamplePhaseRetrieval/data/simdata/testimg.tif";
    %FileTif="C:/research/chisquareeval/data/chivsI9exp_fl/zernike_mode_img/testimg10/roi301zernike2250_13.tif";
else
    FileTif="C:/research/chisquareeval/data/505_515bead60mlambda/"+num2str(imgind)+".tif";
end
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);


for i=1:NumberImages
   im(:,:,i)=imread(FileTif,'Index',i);
end
im=double(im);
im=im./1e6;




parameters.fitmodel='aberrations';

parameters.Nitermax=120;
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


%% Perform fit of uncorrected TFS

parameters.Mx = size(im,1);
parameters.My = size(im,2);
parameters.Mz = size(im,3);
parameters.Ncfg = size(im,4);

parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;

for iii=1:9
    TFSuncorexp=poissrnd(im.*2e3+3);

    [~,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = get_pupil_matrix(parameters);

    [XImage,YImage,~,~] = get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

    % MLE fit
    fprintf('Fitting uncorrected...\n')
    thetainit = initialvalues(TFSuncorexp,XImage,YImage,parameters);
    [thetastore,TFSuncorfit_phase,~,~] = localization(TFSuncorexp,thetainit,parameters);
    thetauncor = squeeze(thetastore(:,:,end));

    for i=4:length(thetauncor)-2
        aberrationsall(i,3)=aberrationsall(i,3)+thetauncor(i);
    end
end

aberrationsall(:,3)=aberrationsall(:,3)/9;
aberrationsp=aberrationsall;
thetauncor_phase=thetauncor;
aberrations=aberrationsall;
parameters.fitmodel='aberrationsamp';
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

[~,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = get_pupil_matrix(parameters);

[XImage,YImage,~,~] = get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

% MLE fit
for iii=1:9
    fprintf('Fitting uncorrected...\n')
    TFSuncorexp=poissrnd(im.*2e3+3);
    thetainit = initialvalues(TFSuncorexp,XImage,YImage,parameters);
    [thetastore,TFSuncorfit,~,numiters] = localization(TFSuncorexp,thetainit,parameters);
    thetauncor = squeeze(thetastore(:,:,end));

    if strcmp(parameters.fitmodel,'aberrationsamp') 
        for i=4:length(aberrations)
            aberrations(i,3)=thetauncor(i);
            aberrations_amp(i,3)=aberrations_amp(i,3)+thetauncor(i);
        end
        for i=1:length(aberrations)
            aberrations(i,4)=thetauncor(length(aberrations)+i);
            aberrations_amp(i,4)=aberrations_amp(i,4)+thetauncor(length(aberrations)+i);
        end
    elseif strcmp(parameters.fitmodel,'aberrations')    
        for i=1:length(aberrations)
            aberrations(i,3)=thetauncor(i+3);
            aberrations_amp(i,3)=aberrations_amp(i,3)+thetauncor(i+3);
        end
    end
end

%save('C:/code/ExamplePhaseRetrieval/data/expdata/vector_fitter_amp','PSF','thetauncor');
%save('C:/code/ExamplePhaseRetrieval/data/expdata/vector_fitter_phase','PSF_phase','thetauncor_phase');
%{
matfile=load('C:/code/ExamplePhaseRetrieval/data/simdata/vector_fitter_amp.mat');
PSFa=matfile.PSF;

matfile=load('C:/code/ExamplePhaseRetrieval/data/simdata/vector_fitter_phase.mat');
PSFp=matfile.PSF_phase;

subplot(3,3,1)
imagesc(PSFp(:,:,5));
subplot(3,3,2)
imagesc(PSFp(:,:,11));
subplot(3,3,3)
imagesc(PSFp(:,:,16));


subplot(3,3,4)
imagesc(PSFa(:,:,5));
subplot(3,3,5)
imagesc(PSFa(:,:,11));
subplot(3,3,6)
imagesc(PSFa(:,:,16));

subplot(3,3,7)
imagesc(im(:,:,5));
subplot(3,3,8)
imagesc(im(:,:,11));
subplot(3,3,9)
imagesc(im(:,:,16));
%}



%{


subplot(3,3,1)
imagesc(TFSuncorexp(:,:,1));
subplot(3,3,2)
imagesc(TFSuncorexp(:,:,10));
subplot(3,3,3)
imagesc(TFSuncorexp(:,:,21));


subplot(3,3,4)
imagesc(TFSuncorfit(:,:,1));
subplot(3,3,5)
imagesc(TFSuncorfit(:,:,10));
subplot(3,3,6)
imagesc(TFSuncorfit(:,:,21));

subplot(3,3,7)
imagesc(TFSuncorfit_phase(:,:,1));
subplot(3,3,8)
imagesc(TFSuncorfit_phase(:,:,10));
subplot(3,3,9)
imagesc(TFSuncorfit_phase(:,:,21));
%}