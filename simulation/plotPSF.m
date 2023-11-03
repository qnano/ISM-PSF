fn="C:/code/jelmer/data/flat_tetrapod1409-100-1.tif";
%fn="C:/code/jelmer/data/0109-2/flat-"+num2str(img_ind)+".tif";
tiff_info = imfinfo(fn); % return tiff structure, one element per image
tiff_stack = imread(fn, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(fn, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

imagesc(tiff_stack(:,:,5));
