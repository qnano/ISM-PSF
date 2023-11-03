clear all
Npupil=201;
NP=11;
N_sam=10;
PSF_em_fitting_arr=zeros([N_sam,NP]);
PSF_im_fitting_arr=zeros([N_sam,NP]);
T_efficiency=zeros([1,NNN]);
%pinhole_radius_at_imaging_plane=640*0.61/1.33; % nm
pinhole_scale=linspace(0.3,2,11);
rms=24.48*1;

for kkk=1:N_sam
for iii=1:NP
disp(iii);
pinhole_radius_at_imaging_plane=640*0.61*pinhole_scale(iii)/1.33; % nm
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

pinhole_radius_at_imaging_plane_pixel=pinhole_radius_at_imaging_plane/parameters.pixelsize;
[xx,yy]=meshgrid(linspace(-n_pixels/2,n_pixels/2,n_pixels),linspace(-n_pixels/2,n_pixels/2,n_pixels));
pinhole=zeros(n_pixels);
pinhole(sqrt(xx.^2+yy.^2)<pinhole_radius_at_imaging_plane_pixel)=1; % pinhole at imaging plane

parameters.zemit =0;

parameters.zrange = [-750, 750];

parameters.Mx = n_pixels ;
parameters.My = n_pixels ;
parameters.Mz = n_z_slices;
parameters.Npupil=Npupil;
parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;
parameters.ztype='medium';
parameters.polarization_excite='circular';  % circular or linear
%aberrations(5,3)=0;
%rms=0;

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


%parameters.Npupil=parameters.Npupil*2;
[OTF_ex] = ...
  get_OTF_matrix(PSF_ex,parameters);

[OTF_pinhole] = ...
  get_OTF_matrix(pinhole,parameters);

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
parameters.Npupil=Npupil;
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


T_efficiency(kkk,iii)=sum(sum(PSF_em.*pinhole))/sum(sum(PSF_em));
%parameters.Npupil=parameters.Npupil*2;

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);
%OTF_im=imaging_scanning_pinhole_OTF(OTF_em,OTF_em,OTF_pinhole,parameters);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im,parameters);
PSF_im=PSF_im./max(max(PSF_im));

[lx,ly]=size(PSF_im);
[XX,YY]=meshgrid(linspace(-parameters.xrange,parameters.xrange,lx),linspace(-parameters.yrange,parameters.yrange,ly));
x_im=fit_2Dgaussian(XX,YY,PSF_im);

x_em=fit_2Dgaussian(XX,YY,PSF_em);

PSF_em_fitting_arr(kkk,iii)=x_em(3);
PSF_im_fitting_arr(kkk,iii)=x_im(3);

end
end

fn="C:\research\spin disk\figure\aberration_and_PSF_width_T_eff\36ml_imaging_scanning.mat";

save(fn,'PSF_im_fitting_arr','PSF_em_fitting_arr','T_efficiency')




