if true
clear all
NNN=11;
PSF_em_fitting_arr=zeros([23,23,NNN]);
PSF_im_fitting_arr=zeros([21,21,NNN]);
T_efficiency=zeros([1,NNN]);


[parameters]=initial_parameter();
pinhole_radius_at_imaging_plane=640*0.61/parameters.NA; % nm
parameters.lambda=640;
parameters.NA=1.3;
n_pixels=parameters.Mx;
pinhole_radius_at_imaging_plane_pixel=pinhole_radius_at_imaging_plane/parameters.pixelsize;
[xx,yy]=meshgrid(linspace(-n_pixels/2,n_pixels/2,n_pixels),linspace(-n_pixels/2,n_pixels/2,n_pixels));
pinhole=zeros(n_pixels);
pinhole(sqrt(xx.^2+yy.^2)<pinhole_radius_at_imaging_plane_pixel)=1; % pinhole at imaging plane

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_ex_pupil_matrix_v3(parameters);

[~,~,FieldMatrix] = ...
  get_ex_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_ex] = get_ex_psfs(FieldMatrix,parameters);

PSF_ex=PSF_ex(:,:,2);

[lx,ly,lz]=size(PSF_ex);
%PSF_ex=PSF_ex(diff:end-diff-1,diff:end-diff-1,3);
PSF_ex=PSF_ex./(sum(sum(PSF_ex)));
%parameters.Npupil=parameters.Npupil*2;
[OTF_ex] = ...
  get_OTF_matrix(PSF_ex,parameters);

[OTF_pinhole] = ...
  get_OTF_matrix(pinhole,parameters);


[parameters]=initial_parameter();
parameters.lambda=690;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
[lx,ly,lz]=size(PSF_em);
diff=(lx-n_pixels)/2;
PSF_em=PSF_em(diff:end-diff-1,diff:end-diff-1,3);

PSF_em=PSF_em./(sum(sum(PSF_em)));

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im,parameters);
PSF_im=PSF_im./max(max(PSF_im));


fn="D:\spin disk\figure new\figure1\OTFandPSF_vectorial_stokesshift.mat";
save(fn,'PSF_em','OTF_em','OTF_im','PSF_im','PSF_ex','OTF_ex')

end

if true
%%%%%%%%%%
%scalar+no Stokes shift
%%%%%%%%%%

clear all


rms=24.48*0;
[parameters]=initial_parameter();
pinhole_radius_at_imaging_plane=640*0.61/parameters.NA; % nm
parameters.lambda=640;
parameters.NA=1.3;
n_pixels=parameters.Mx;

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
diff=(lx-n_pixels)/2;
%PSF_ex=PSF_ex(diff:end-diff-1,diff:end-diff-1,3);

PSF_ex=PSF_ex./(sum(sum(PSF_ex)));


%parameters.Npupil=parameters.Npupil*2;
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


%T_efficiency(iii)=sum(sum(PSF_em.*pinhole))/sum(sum(PSF_em));
%parameters.Npupil=parameters.Npupil*2;

[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

%OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);
OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im,parameters);
PSF_im=PSF_im./max(max(PSF_im));

fn="D:\spin disk\figure new\figure1\OTFandPSF_scalar.mat";
save(fn,'PSF_em','OTF_em','OTF_im','PSF_im','PSF_ex','OTF_ex')

end

if true
%%%%%%%%%%
%vectorial+no Stokes shift
%%%%%%%%%%

clear all
NNN=11;


[parameters]=initial_parameter();
pinhole_radius_at_imaging_plane=640*0.61/parameters.NA; % nm
parameters.lambda=640;
parameters.NA=1.3;
n_pixels=parameters.Mx;

pinhole_radius_at_imaging_plane_pixel=pinhole_radius_at_imaging_plane/parameters.pixelsize;
[xx,yy]=meshgrid(linspace(-n_pixels/2,n_pixels/2,n_pixels),linspace(-n_pixels/2,n_pixels/2,n_pixels));
pinhole=zeros(n_pixels);
pinhole(sqrt(xx.^2+yy.^2)<pinhole_radius_at_imaging_plane_pixel)=1; % pinhole at imaging plane

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_ex_pupil_matrix_v3(parameters);

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

[parameters]=initial_parameter();
parameters.lambda=640;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
[lx,ly,lz]=size(PSF_em);
diff=(lx-n_pixels)/2;
PSF_em=PSF_em(diff:end-diff-1,diff:end-diff-1,3);

PSF_em=PSF_em./(sum(sum(PSF_em)));


[OTF_em] = ...
  get_OTF_matrix(PSF_em,parameters);

OTF_pinhole=OTF_pinhole./max(max(OTF_pinhole));

OTF_im=imaging_scanning_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters);

[PSF_im] = get_PSFfromOTF_matrix(OTF_im,parameters);
PSF_im=PSF_im./max(max(PSF_im));

fn="D:\spin disk\figure new\figure1\OTFandPSF_vectorial_no_stokesshift.mat";
save(fn,'PSF_em','OTF_em','OTF_im','PSF_im','PSF_ex','OTF_ex')
end
