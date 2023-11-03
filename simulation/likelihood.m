function [logL,gradlogL,HessianlogL] = likelihood(M,mu,dmudtheta,varfit)
% returns the log-likelihood, as well as first and second order
% derivatives w.r.t. the parameters for a noisy image M, measured in
% number of photons per pixel, and Poisson-rate mu with first order
% derivatives dmudtheta. The log-likelihood is Poisson+readout-noise based.
% The script can work for a single ROI data/model: size(mu)=[Mx,My] or for
% multiple ROI data/model: size(mu)=[Mx,My,Mz].
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%

ndims = length(size(mu));
numparams = size(dmudtheta,ndims+1);

% calculation of weight factors
keps = 1e3*eps;
mupos = double(mu>0).*mu + double(mu<0)*keps;
weight = (M-mupos)./(mupos+varfit);
dweight = (M+varfit)./(mupos+varfit).^2;

if ndims==2
% log-likelihood merit function
  logL = sum(sum( (M+varfit).*log(mupos+varfit)-(mupos+varfit) ));

% gradient vector and Hessian matrix
  gradlogL = zeros(1,numparams);
  HessianlogL = zeros(numparams,numparams);
  for ii = 1:numparams
    gradlogL(ii) = sum(sum( weight.*dmudtheta(:,:,ii) ));
    for jj = 1:numparams
      HessianlogL(ii,jj) = sum(sum(-dweight.*dmudtheta(:,:,ii).*dmudtheta(:,:,jj)));
    end
  end
end

if ndims==3
% log-likelihood merit function
  logL = sum(sum(sum( (M+varfit).*log(mupos+varfit)-(mupos+varfit) )));

% gradient vector and Hessian matrix
  gradlogL = zeros(1,numparams);
  HessianlogL = zeros(numparams,numparams);
  for ii = 1:numparams
    gradlogL(ii) = sum(sum(sum( weight.*dmudtheta(:,:,:,ii) )));
    for jj = 1:numparams
      HessianlogL(ii,jj) = sum(sum(sum(-dweight.*dmudtheta(:,:,:,ii).*dmudtheta(:,:,:,jj))));
    end
  end
end

return