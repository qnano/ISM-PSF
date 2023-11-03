
NNN=11;
PSF_em_fitting_arr=zeros([23,23,NNN]);
PSF_im_fitting_arr=zeros([21,21,NNN]);

for iii=1:NNN
load('C:/research/issue1/simulator/p.mat');
count=1;
ma = 6;
for i=1:ma
    mind=-i:2:i;
    for j=1:length(mind)
        aberrations(count,:)=[i,mind(j),0,0];
        count=count+1;
    end
end
aberrations(3,2)=0;
aberrations(4,2)=-2;

n_pixels =512;

n_z_slices=3;
parameters.lambda=640;
parameters.pixelsize=5;
parameters.zemit =0;

parameters.zrange = [-750, 750];

parameters.Mx = n_pixels ;
parameters.My = n_pixels ;
parameters.Mz = n_z_slices;
parameters.Npupil=512;
parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;
parameters.ztype='medium';
parameters.polarization_excite='circular';  % circular or linear
%aberrations(5,3)=0;
rms=24.48*5;
aberr=rand(24,1)-0.5;
aberr=aberr./(sqrt(sum(aberr.^2)));
aberr=aberr.*rms;
aberrations(4:end,3)=aberr;
parameters.aberrations=aberrations;
parameters.numparams=5+length(aberrations);


parameters.xemit = 0;
parameters.yemit =0;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_ex_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_ex_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_ex] = get_ex_psfs(FieldMatrix,parameters);

PSF_ex=PSF_ex(:,:,2);

[lx,ly,lz]=size(PSF_ex);
diff=(lx-n_pixels)/2;
%PSF_ex=PSF_ex(diff:end-diff-1,diff:end-diff-1,3);

PSF_ex=PSF_ex./(sum(sum(PSF_ex)));


parameters.Npupil=parameters.Npupil*2;
[OTF_ex] = ...
  get_OTF_matrix(PSF_ex,parameters);


load('C:/research/issue1/simulator/p.mat');
count=1;
ma = 6;
for i=1:ma
    mind=-i:2:i;
    for j=1:length(mind)
        aberrations(count,:)=[i,mind(j),0,0];
        count=count+1;
    end
end
aberrations(3,2)=0;
aberrations(4,2)=-2;



n_z_slices=3;

parameters.lambda=690;


parameters.zemit =0;

parameters.zrange = [-750, 750];

parameters.Mx = n_pixels ;
parameters.My = n_pixels ;
parameters.Mz = n_z_slices;
parameters.Npupil=512;
parameters.pixelsize=5;

parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;
parameters.ztype='medium';
aberrations(4:end,3)=aberr;
parameters.aberrations=aberrations;
parameters.numparams=5+length(aberrations);


parameters.xemit = 0;
parameters.yemit =0;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
[lx,ly,lz]=size(PSF_em);
diff=(lx-n_pixels)/2;
PSF_em=PSF_em(diff:end-diff-1,diff:end-diff-1,3);

PSF_em=PSF_em./(sum(sum(PSF_em)));

parameters.Npupil=parameters.Npupil*2;

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_im=imaging_scanning_OTF(OTF_ex,OTF_em,parameters);

l=sqrt(length(OTF_im));

OTF_im_sqrt=reshape(OTF_im,[l,l]);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im_sqrt,parameters);
PSF_im=PSF_im./max(max(PSF_im));


[lx,ly]=size(PSF_im);
%[XX,YY]=meshgrid(linspace(-parameters.xrange,parameters.xrange,lx),linspace(-parameters.yrange,parameters.yrange,ly));
%x_im=fit_2Dgaussian(XX,YY,PSF_im);

%x_em=fit_2Dgaussian(XX,YY,PSF_em);

%disp(x_em(3)/x_im(3))



parameters.Mx = 21 ;
parameters.My = 21 ;

parameters.pixelsize=65;

parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;

[PSF_im_fitting] = get_PSFfromOTF_matrix(OTF_im_sqrt,parameters);

parameters.lambda=690;
[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em_fitting] = get_psfs(FieldMatrix,parameters);
PSF_em_fitting=PSF_em_fitting(:,:,3);

PSF_em_fitting_arr(:,:,iii)=PSF_em_fitting;
PSF_im_fitting_arr(:,:,iii)=PSF_im_fitting;

end

fn="C:\research\spin disk\simulation\PSF\PSFs_aberration_180ml.mat";
save(fn,'PSF_im_fitting_arr','PSF_em_fitting_arr')
if false
figure()

PupilSize = 4.0;
    
DxyPupil = 2*PupilSize/parameters.Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
plot(XYPupil,abs(OTF_em(256,:))./max(abs(OTF_em(256,:))));
hold on
plot(XYPupil,abs(OTF_im_sqrt(256,:))./max(abs(OTF_im_sqrt(256,:))));

legend('WF OTF','IM OTF')
end
if false
PSF_im=PSF_im(256-100:256+100,256-100:256+100);
show_em=PSF_em(256-100:256+100,256-100:256+100);
show_em=show_em./max(max(show_em));
all=horzcat(PSF_im,show_em);

figure()
imagesc(all)
end



