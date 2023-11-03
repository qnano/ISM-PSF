clear all
NNN=11;
PSF_em_fitting_arr=zeros([23,23,NNN]);
PSF_im_fitting_arr=zeros([21,21,NNN]);
T_efficiency=zeros([1,NNN]);


[parameters]=initial_parameter();
parameters.polarization_excite='linear';
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
save("xlinear_expsf.mat",'PSF_ex');


clear all
NNN=11;
PSF_em_fitting_arr=zeros([23,23,NNN]);
PSF_im_fitting_arr=zeros([21,21,NNN]);
T_efficiency=zeros([1,NNN]);


[parameters]=initial_parameter();
parameters.polarization_excite='ylinear';
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
save("ylinear_expsf.mat",'PSF_ex');

clear all
NNN=11;
PSF_em_fitting_arr=zeros([23,23,NNN]);
PSF_im_fitting_arr=zeros([21,21,NNN]);
T_efficiency=zeros([1,NNN]);


[parameters]=initial_parameter();
parameters.polarization_excite='circular';
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
save("circular_expsf.mat",'PSF_ex');


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
save("scalar_expsf.mat",'PSF_ex');



clear all

[parameters]=initial_parameter();
parameters.lambda=690;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
save("vectorial_empsf.mat",'PSF_em');


clear all

[parameters]=initial_parameter();
parameters.lambda=690;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
PSF_em=PSF_em(:,:,3);
PSF_em=PSF_em./(sum(sum(PSF_em)));

save("vectorial_empsf.mat",'PSF_em');

clear all
[parameters]=initial_parameter();
parameters.lambda=640;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_pupil_matrix_scalar(parameters);

[~,~,FieldMatrix] = ...
  get_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_em] = get_psfs(FieldMatrix,parameters);
PSF_em=PSF_em(:,:,3);
PSF_em=PSF_em./(sum(sum(PSF_em)));

save("scalar_empsf.mat",'PSF_em');
