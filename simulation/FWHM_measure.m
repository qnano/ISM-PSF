function FWHM=FWHM_measure(psf)
coor=1:1:512;
coorfine=1:0.2:512;

x=psf(256,:);
xfine = interp1(coor,x,coorfine);
ma=max(xfine);

for i=1:length(xfine)
if xfine(i)>0.5*ma
    x1=i;
    break;
end
end


for i=1:length(xfine)
if xfine(length(xfine)-i+1)>0.5*ma
    x2=length(xfine)-i+1;
    break;
end
end

FWHM=x2-x1;

end