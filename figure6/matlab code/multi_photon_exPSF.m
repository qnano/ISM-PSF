clear all

[parameters]=initial_parameter();
n_pixels=parameters.Mx;

pinhole_radius_at_imaging_plane=640*0.61/parameters.NA; % nm
pinhole_radius_at_imaging_plane_pixel=pinhole_radius_at_imaging_plane/parameters.pixelsize;
[xx,yy]=meshgrid(linspace(-n_pixels/2,n_pixels/2,n_pixels),linspace(-n_pixels/2,n_pixels/2,n_pixels));
pinhole=zeros(n_pixels);
pinhole(sqrt(xx.^2+yy.^2)<pinhole_radius_at_imaging_plane_pixel)=1; % pinhole at imaging plane

PSF_ex_arr=zeros([3,512,512]);

for i=1:3
parameters.lambda=640.0*i;

[tt,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
         get_ex_pupil_matrix_v3(parameters);

[~,~,FieldMatrix] = ...
  get_ex_field_matrix(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);

[PSF_ex] = get_ex_psfs(FieldMatrix,parameters);
PSF_ex_arr(i,:,:)=abs(PSF_ex(:,:,2));
end



fn="D:\spin disk\figure new\figure6\multi_photon_exPSF.mat";
save(fn,'PSF_ex_arr')