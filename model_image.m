    rng('shuffle');
    sigl=1.15*65;%nm
    sigh=1.9*65;%nm
    th=800;% the range of commet tail decay length, max (unit:nm)
    tl=80;% the range of commet tail decay length, min (unit:nm)
    Num=60;
    x=1:1:6*Num*5; % 13nm per pixel
    A=rand(2000,1);
    B=rand(2000,1);
    %% prepare initial
    jsq=0;
    t=((th-tl)*B(jsq+500)+tl)/13;% 13nm pixel
    sigma=((sigh-sigl)*A(jsq+400)+sigl)/13;
    mu_est=50*sigma;
    gauss_psf=norm_density(x, mu_est, sigma);
    exp_decay= zeros(length(x),1); 
    exp_decay(1:round(50*sigma))=0;
    exp_decay(round(50*sigma):end)=exp(-(x(round(50*sigma):end)-round(50*sigma))/t);
    exp_gauss=conv(exp_decay,gauss_psf);
    mask=ones(5,1);%13nm/pixel to 65nm/pixel
    exp_gauss_c=conv(exp_gauss,mask,'valid');
    eg_c=exp_gauss_c(1:5:end);
    egc=eg_c/max(eg_c);
    eg_data=egc*150;% 150 is a tunable parameter, determined based on the intensity range of the signal from your experiment
    alpha=0.75; %photon number

    %% add photon shot noise
avg_shot=poissrnd(eg_data,1,length(eg_data));
%% alpha conversion factor
avg_read=uint16((avg_shot*alpha)/0.53+100+normrnd(0,2,[1,length(eg_data)]));% add readout noise
% avg_read=uint16((eg_data*alpha)/0.53);
[maxa,maxind]=max(avg_read);
data=avg_read(maxind-20:maxind+50);
profile_meanI=double(data);

[pmax,pI]=max(profile_meanI);

save data.mat