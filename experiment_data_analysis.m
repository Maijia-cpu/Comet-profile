load('trackdata.mat')
% trackdata.mat contains profile_meanI (7 pixel width sum of intensity profile against comet moving direction, background is subtracted from profile_meanI)
alpha=0.75;%quantum efficiency of the SCMOS camera at wavelength 488nm
profile_meanI=flip(profile_meanI);% change to along comet moving direction
avg_read=(profile_meanI)*0.53/alpha;% convert ADU to photon,0.53 is the CMOS sensitivity (unit e/ADU)
maxI=max(avg_read);
[pmax,pI]=max(profile_meanI);
%% run the filtting algorithm six times
parfor i=1:6
[param(:,:,i), profiledata(:,:,i), meanbg(:,i), bgvalue(:,:,i),lambda_value(:,i),sigmavalue(:,i),coeff_amp(:,i),center_value(:,i)] = fit_conv_new(profile_meanI,pI,maxI);
end


for j=1:6
psfsigma(j)=mean(param(end-200:end,1,j))*13;% unit nm standard devaition of point spread function 
tail(j)=mean(param(end-200:end,2,j))*13;% decay length
end