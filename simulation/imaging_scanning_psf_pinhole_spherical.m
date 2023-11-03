clear all
Npupil=201;
NNN=11;
PSF_em_fitting_arr=zeros([33,33,NNN]);
PSF_im_fitting_arr=zeros([31,31,NNN]);

pinhole_radius_at_imaging_plane=640*0.61/1.33; % nm
fn="C:/research/spin disk/figure/index_mismatch_crlb/index_mismatch/glycerol-water/sph_arr2.mat";

mat=load(fn);
sph_arr=mat.sph_arr;
T_efficiency=zeros([1,NNN]);

for iii=1:NNN
disp(iii);

load('C:/research/issue1/simulator/p.mat');
aberrations(1,:)=[4,0,sph_arr(iii,1),0];
aberrations(2,:)=[6,0,sph_arr(iii,2),0];
aberrations(3,:)=[8,0,sph_arr(iii,3),0];
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
T_efficiency(iii)=sum(sum(PSF_em.*pinhole))/sum(sum(PSF_em));
%parameters.Npupil=parameters.Npupil*2;

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);

%[XX,YY]=meshgrid(linspace(-parameters.xrange,parameters.xrange,lx),linspace(-parameters.yrange,parameters.yrange,ly));
%x_im=fit_2Dgaussian(XX,YY,PSF_im);

%x_em=fit_2Dgaussian(XX,YY,PSF_em);

%disp(x_em(3)/x_im(3))



parameters.Mx = 31 ;
parameters.My = 31 ;

parameters.pixelsize=65;

parameters.xrange = parameters.pixelsize*parameters.Mx/2;
parameters.yrange = parameters.pixelsize*parameters.My/2;

[PSF_im_fitting] = get_PSFfromOTF_matrix(OTF_im,parameters);

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

fn="C:/research/spin disk/figure/index_mismatch_crlb/index_mismatch/glycerol-water/PSFs_spherical_aberration_pinhole1.mat";

save(fn,'PSF_im_fitting_arr','PSF_em_fitting_arr','T_efficiency')

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



