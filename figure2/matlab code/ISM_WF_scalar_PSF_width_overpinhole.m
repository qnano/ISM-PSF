clear all

NNN=31;
PSF_em_fitting_arr=zeros([1,NNN]);
PSF_im_fitting_arr=zeros([1,NNN]);
T_efficiency=zeros([1,NNN]);
pinhole_scale=linspace(0.3,2,NNN);


for iii=1:NNN
disp(iii);

[parameters]=initial_parameter();
n_pixels=parameters.Mx;

pinhole_radius_at_imaging_plane=640*0.61*pinhole_scale(iii)/parameters.NA; % nm
parameters.lambda=640;

pinhole_radius_at_imaging_plane_pixel=pinhole_radius_at_imaging_plane/parameters.pixelsize;
[xx,yy]=meshgrid(linspace(-n_pixels/2,n_pixels/2,n_pixels),linspace(-n_pixels/2,n_pixels/2,n_pixels));
pinhole=zeros(n_pixels);
pinhole(sqrt(xx.^2+yy.^2)<pinhole_radius_at_imaging_plane_pixel)=1; % pinhole at imaging plane

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_ex_pupil_matrix_scalar(parameters);

[~,~,FieldMatrix] = ...
  get_ex_field_matrix_scalar(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_ex] = get_ex_psfs_scalar(FieldMatrix,parameters);

PSF_ex=PSF_ex(:,:,2);

[lx,ly,lz]=size(PSF_ex);

PSF_ex=PSF_ex./(sum(sum(PSF_ex)));

[OTF_ex] = ...
  get_OTF_matrix(PSF_ex,parameters);

[OTF_pinhole] = ...
  get_OTF_matrix(pinhole,parameters);

[parameters]=initial_parameter();
parameters.lambda=640;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix_scalar(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
[lx,ly,lz]=size(PSF_em);
diff=(lx-n_pixels)/2;
PSF_em=PSF_em(diff:end-diff-1,diff:end-diff-1,3);

PSF_em=PSF_em./(sum(sum(PSF_em)));


T_efficiency(iii)=sum(sum(PSF_em.*pinhole))/sum(sum(PSF_em));

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

%OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);
OTF_im=imaging_scanning_pinhole_OTF(OTF_em,OTF_em,OTF_pinhole,parameters);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im,parameters);
PSF_im=PSF_im./max(max(PSF_im));


[lx,ly]=size(PSF_im);
[XX,YY]=meshgrid(linspace(-parameters.xrange,parameters.xrange,lx),linspace(-parameters.yrange,parameters.yrange,ly));
x_im=fit_2Dgaussian(XX,YY,PSF_im);

x_em=fit_2Dgaussian(XX,YY,PSF_em);

PSF_em_fitting_arr(iii)=x_em(3);
PSF_im_fitting_arr(iii)=x_im(3);


end

fn="D:\spin disk\figure new\figure2\new\PSFs_width_oversizeof pinhole.mat";

save(fn,'PSF_im_fitting_arr','PSF_em_fitting_arr','pinhole_scale')




