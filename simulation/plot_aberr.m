figure()
x=1:length(aberrations_amp);
name={};
for i=1:length(aberrations_amp)
    
    name=cat(2,name,{num2str(aberrations_amp(i,1))+","+num2str(aberrations_amp(i,2))});
end
y = aberrations_amp(:,3);
plot(x,y,"--*")
xticks(x)
xticklabels(name)
xtickangle(45)
xlabel("(n,m)");
ylabel("aberration (m\lambda)");
figure()
Npupil=65;
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

orders = aberrations_amp(:,1:2);
zernikecoefs = squeeze(parameters.aberrations(:,3));
normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
zernikecoefs = normfac.*zernikecoefs;
% allzernikes is normalized
allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
Waberration = zeros(size(XPupil));
for j = 1:numel(zernikecoefs)
  Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));  
end
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
Waberration=Waberration.*ApertureMask;
colormap jet
imagesc(Waberration);
title("pupil phase");
colorbar()