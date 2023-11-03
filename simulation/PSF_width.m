fn="C:/research/spin disk/figure/index_mismatch_crlb/index_mismatch/oil-water/PSF_width.mat";
mat=load(fn);
PSF_em_image=mat.PSF_em_image;
PSF_im_image=mat.PSF_im_image;

xrange=512*5/2;
yrange=512*5/2;
[lx,ly]=size(PSF_im_image(:,:,1));
[XX,YY]=meshgrid(linspace(-xrange,xrange,lx),linspace(-yrange,yrange,ly));
for ii=1:11
x_im=fit_2Dgaussian(XX,YY,PSF_im_image(:,:,ii));

x_em=fit_2Dgaussian(XX,YY,PSF_em_image(:,:,ii));

PSF_em_fitting_arr(ii)=x_em(3);
PSF_im_fitting_arr(ii)=x_im(3);
end