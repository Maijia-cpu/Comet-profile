%% fit  exponential decays with convolution
function [param, profiledata,mean_bg,bgvalue,lambda_value, sigmavalue,coeff_amp,center_value] = fit_conv_new(profile_meanI,pIo,maxI)
Num=length(profile_meanI)-1;
xm=(0:1:Num)*0.065;% 0.065micron per pixel 
pmax=max(profile_meanI);

%% 1. fit two exponential decays first to extract parameter fitting ranges for the following fittings
lambda=zeros(2,10);
sum_all=zeros(10,1);
% Uses fitnlm() to fit a non-linear model (an exponential decay curve) through noisy data.
% Requires the Statistics and Machine Learning Toolbox, which is where fitnlm() is contained.
% Initialization steps.
    for ip=1:14
    clear coefficients
    % counter motion direction
    % fit the first half
    Xc = xm(1:round(pIo-(ip-7.5)*0.5));
    Yc = profile_meanI(1:round(pIo-(ip-7.5)*0.5)); 
%     plot(Xc, Yc, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
%     grid on;
    % Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
    tbl = table(Xc', Yc');
    % Define the model as Y = a + exp(-b*x)
    % Note how this "x" of modelfun is related to big X and big Y.
    modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*(-x(:, 1)+(pIo-(ip-7.5)*0.5-1)*0.065));  
    beta0 = [profile_meanI(1), pmax-profile_meanI(1), 1]; % Guess values to start with.  Just make your best guess.
    % Now the next line is where the actual model computation is done.
    try
    mdl = fitnlm(tbl, modelfun, beta0);
    coefficients = mdl.Coefficients{:, 'Estimate'};
    catch ME
           switch ME.identifier
               case 'stats:nlinfit:NonFiniteFunOutput'
                 lambda(1,ip)=Inf;  
                 sumdc=Inf;
           end
    end
    % Extract the coefficient values from the the model object.
    % The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
    Ac=exist('coefficients');
    if Ac>0
    % Create smoothed/regressed data using the model:
    yFitted = coefficients(1) + coefficients(2) * exp(-coefficients(3)*(-Xc+(pIo-(ip-7.5)*0.5-1)*0.065));
    % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
%     hold on;
%     plot(Xc, yFitted, 'r-', 'LineWidth', 2);
%     grid on;
    sumdc=sum((yFitted-Yc).^2);
    lambda(1,ip)=1/coefficients(3);
    end
    clear coefficients
    % along motion direction
    % fit the second half
    Xc = xm(round(pIo-(ip-7.5)*0.5)+1:end);
    Yc = profile_meanI(round(pIo-(ip-7.5)*0.5)+1:end); 
%     plot(Xc, Yc, 'b*', 'LineWidth', 2, 'MarkerSize', 15);
%     grid on;
    tbl = table(Xc', Yc');
    modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*(x(:, 1)+(-pIo+(ip-7.5)*0.5)*0.065));  
    beta0 = [profile_meanI(1), pmax-profile_meanI(1), 1]; % Guess values to start with.  Just make your best guess.
    % Now the next line is where the actual model computation is done.
   try
    mdl = fitnlm(tbl, modelfun, beta0);
    coefficients = mdl.Coefficients{:, 'Estimate'};
    catch ME
           switch ME.identifier
               case 'stats:nlinfit:NonFiniteFunOutput'
                 lambda(1,ip)=Inf;  
                 sumdp=Inf;
           end
   end
    AB=exist('coefficients');
    if AB>0
    % Create smoothed/regressed data using the model:
    yFitted = coefficients(1) + coefficients(2) * exp(-coefficients(3)*(Xc+(-pIo+(ip-7.5)*0.5)*0.065));
    % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
%     hold on;
%     plot(Xc, yFitted, 'r-', 'LineWidth', 2);
%     grid on;
    sumdp=sum((yFitted-Yc).^2);
    lambda(2,ip)=1/coefficients(3);  
    end
    sum_all(ip)=sumdp+sumdc;% sum of squared difference
    end
    
    [smin,minInd]=min(sum_all);
    lambda_l=lambda(:,minInd);% the fitted decay lengths of two exponential decays
    bg1=profile_meanI(1:round(pIo-lambda_l(1)*1000/65*2)-1);%background of one side 
    bg2=profile_meanI(round(pIo+lambda_l(2)*1000/65*2)+1:end);% background of the other side
    bg=[bg1,bg2];
    mean_bg=mean(bg);% background value 
%      if lambda_l(1)>lambda_l(2)% reverse data
%          profile_meanI=flip(profile_meanI);
%          fliprecord=1;
%      else
%          fliprecord=0;
%      end  
     profilemeanI=profile_meanI-mean_bg;% subtract background
     
%% 2. fit left side with gaussian and right side with exponential decay to extract parameter ranges for analytical fittings
       pI_data=round(pIo-(minInd-7.5)*0.5);% peak position determined from two exponential decay fits
%        if fliprecord==1
%            pI=length(profile_meanI)+1-pI_data;
%        else
           pI=pI_data;
%        end
      
    for ig=1:11
       gaus_profile=profilemeanI(1:pI-(ig-6));
       Ng=length(gaus_profile);
       xg=(0:1:Ng-1)*0.065;
       clear mycoeff
        try
       f = fit(xg',gaus_profile','gauss1');
       mycoeff_f=coeffvalues(f);
       gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
       startPoints = [mycoeff_f(1) mycoeff_f(2) mycoeff_f(3) 0];
       [f1,gof] = fit(xg',gaus_profile',gaussEqn,'Start', startPoints);
       mycoeff=coeffvalues(f1);
       sigma_gaus(ig)=mycoeff(3)/sqrt(2)*1000;
       ygaus=mycoeff(1)*exp(-((xg-mycoeff(2))/mycoeff(3)).^2)+mycoeff(4);
         catch ME
           switch ME.identifier
               case 'curvefit:fit:infComputed'
                lambda_exp(ig)=Inf;
                sumdif(ig)=Inf; 
               case 'curvefit:fit:nanComputed'
                lambda_exp(ig)=Inf;
                sumdif(ig)=Inf;  
               otherwise
                rethrow(ME)
           end
        end 
         A=exist('mycoeff');
       if (A>0)
         center_gaus=mycoeff(2)/0.065;
        if (mycoeff(1)>0 && mycoeff(2)>0)&& (center_gaus<Num)
%        figure,plot(xg,gaus_profile)
%        hold on
%        plot(xg,ygaus)
           if length(ygaus)>round(center_gaus)+1 
           ydif_gaus=sum((gaus_profile(1:round(center_gaus)+1)-ygaus(1:round(center_gaus)+1)).^2);
           else
           ydif_gaus1=sum((gaus_profile-ygaus).^2); 
           xd=(0:1:round(center_gaus))*0.065;
           ygaus2=mycoeff(1)*exp(-((xd-mycoeff(2))/mycoeff(3)).^2)+mycoeff(4);
           ydif_gaus2=sum((ygaus2(length(ygaus)+1:round(center_gaus)+1)-profilemeanI(length(ygaus)+1:round(center_gaus)+1)).^2);
           ydif_gaus=ydif_gaus1+ydif_gaus2;
           end
%% fit exponential decay based on center from gaussian fitting
       clear coefficients
        Yc = profilemeanI(round(center_gaus)+2:end); 
        Ngc=length(Yc);
        Xc=(0:1:Ngc-1)*0.065;
%         figure,plot(Xc, Yc, 'b*', 'LineWidth', 2, 'MarkerSize', 15);
%         grid on;
        tbl = table(Xc', Yc');
        modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*(x(:, 1)-b(4)));  
        beta0 = [profilemeanI(end), max(profilemeanI)-profilemeanI(end), 1,0]; % Guess values to start with.  Just make your best guess.
        % Now the next line is where the actual model computation is done.
        mdl = fitnlm(tbl, modelfun, beta0);
        coefficients = mdl.Coefficients{:, 'Estimate'};
        % Create smoothed/regressed data using the model:
        yFitted = coefficients(1) + coefficients(2) * exp(-coefficients(3)*(Xc-coefficients(4)));
        % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
%         hold on;
%         plot(Xc, yFitted, 'r-', 'LineWidth', 2);
%         grid on;
        ydif_exp=sum((yFitted-Yc).^2);
        
        centervalue(ig)=center_gaus;
        lambda_exp(ig)=1/coefficients(3);
        sumdif(ig)=ydif_exp+ydif_gaus;
           else
        centervalue(ig)=Inf;
        lambda_exp(ig)=Inf;
        sumdif(ig)=Inf;
         end
       else
        centervalue(ig)=Inf;
        lambda_exp(ig)=Inf;
        sumdif(ig)=Inf;
       end
       clear coefficients
    end
   
       [min_s,minInd_s]=min(sumdif);
       lambda_ep=lambda_exp(minInd_s);%decay length 
       sigma_gaus_ep=sigma_gaus(minInd_s);% standard deviation of gaussian distribution

%% 3. fit with analytical solution using parameters from fitting step #2
    profilemeanI=profilemeanI/max(profilemeanI);% normalized
    xa=(1:1:Num+1)*0.065;
    startcoeffa  = max(profilemeanI)/2;
    startlambda  = 1/lambda_ep;
    startconst   = (pI-(minInd_s-6)-1)*0.065;
    
    if sigma_gaus_ep-100>1.15*65
      sigl=sigma_gaus_ep-100; 
    else
      sigl=1.15*65;%nm determined from tetra speck beads measurement mean-2SD
    end
    startsigma  = sigl/1000;
    startcoeffb  = 0;  
%         figure,plot(Xc, Yc, 'b*', 'LineWidth', 2, 'MarkerSize', 15);
%         grid on;
%% fit analyitical solution based on parameter ranges from fitting 
        tbl = table(xa', profilemeanI');
        modelfun = @(b,x) b(1) * exp(-b(2)*((x-b(3))-b(4)^2*b(2)/2)).* (1- erf((b(4)^2*b(2)-(x-b(3)))/sqrt(2)/b(4)))+b(5);  
        beta0 = [startcoeffa startlambda startconst startsigma...
        startcoeffb]; % Guess values to start with.  Just make your best guess.
        % Now the next line is where the actual model computation is done.
        mdl = fitnlm(tbl, modelfun, beta0);
        coefficients_ana = mdl.Coefficients{:, 'Estimate'};
        % Create smoothed/regressed data using the model:
%         yFitted = coefficients_ana(1) * exp(-coefficients_ana(2)*((x-coefficients_ana(3))-coefficients_ana(4)^2*coefficients_ana(2)/2)).* (1- erf((coefficients_ana(4)^2*coefficients_ana(2)-(x-coefficients_ana(3)))/sqrt(2)/coefficients_ana(4)))+coefficients_ana(5);  
        lambda_value=1/coefficients_ana(2)*1000;% commet decay length
        center_value=round(coefficients_ana(3)/0.065);% center position of the commet
        bgvalue_record=coefficients_ana(5);% background value
        sigmavalue=coefficients_ana(4)*1000;%standard deviation of point spread function
        coeff_amp=coefficients_ana(1);
        
%% MC simulation based on fitted parameters from analytical fitting result
% MC simulation is optional if you are interested in decay length of the
% comet profile
% our experience is that MC simulation can give better prediction of point
% spread function width
 % random number
    jsqn=0;
    while jsqn<1500
    number=floor(2*rand);
    jsqn=jsqn+1;
    randnum(jsqn)=number;
    end
    rng('default');
    Ar=rand(2000,1);
    B=rand(2000,1);
    C=rand(2000,1);
    D=rand(2000,1);
    
    if sigmavalue-100>1.9*65
      sigl=1.15*65; %nm determined from tetra speck beads measurement mean-2SD
    elseif sigmavalue-100<1.9*65
      sigl=sigmavalue-100; 
    end
    
    if sigmavalue+50<1.9*65
      sigh=sigmavalue+50;
    else
      sigh=1.85*65;%nm determined from tetra speck beads measurement mean+2SD
    end
    
    th=(lambda_value+60);
    if (lambda_value-60)>20
    tl=(lambda_value-60);
    else
    tl=20;
    end
    
    x=1:1:6*Num*5; % 13nm per pixel

    % locate peak of data and peak of exp_gauss
    l1=round(2*lambda_l(1)*1000/13)+1;% from two exponential decay fittings
    l2=round(2*lambda_l(2)*1000/13)+1;
     
    %% prepare initial
    jsq=0;
    t=((th-tl)*B(jsq+500)+tl)/13;% 13nm pixel
    sigma=((sigh-sigl)*Ar(jsq+400)+sigl)/13;% 13nm pixel
    mu_est=50*sigma;
    gauss_psf=norm_density(x, mu_est, sigma);
    exp_decay= zeros(length(x),1); 
    exp_decay(1:round(50*sigma))=0;
    exp_decay(round(50*sigma):end)=exp(-(x(round(50*sigma):end)-round(50*sigma))/t);
    exp_gauss=conv(exp_decay,gauss_psf);
    [h,profile,bgdata]=sqsum(exp_gauss,profilemeanI,l1,l2,center_value,bgvalue_record,maxI);
    %% set temperature
    for jt=1:1:1500
    if jt<200
    T(jt)=(0.2-jt*0.1/200)*0.2;  %%%tunable
    elseif (jt>199)
    T(jt)=(0.1-0.06/1600*(jt-200))*0.2;
    end
    end
%% main part
      while jsq<1000
        num=randnum(jsq+1);% random generate 0 or 1 or 2; 
        if num==0
            jsq=jsq+1;
            sigma1=((sigh-sigl)*Ar(jsq+400)+sigl)/13;%13nm/pixel
            mu_est=50*sigma1;
            gauss_psf=norm_density(x, mu_est, sigma1);
            exp_decay= zeros(length(x),1); 
            exp_decay(1:round(50*sigma1))=0;
            exp_decay(round(50*sigma1):end)=exp(-(x(round(50*sigma1):end)-round(50*sigma1))/t);
            exp_gauss=conv(exp_decay,gauss_psf);
            [h1,profile1,bgdata1]=sqsum(exp_gauss,profilemeanI,l1,l2,center_value,bgvalue_record,maxI);
            Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
            % param(:,1) standard deviation of point spread function
            % param(:,2) decay length of commet tail
            % param(:,3) the difference of simulated image and experimental
            % image
            if Z<1
                num1=B(500+jsq);
                if Z>num1
                    sigma=sigma1;
                    h=h1;
                    param(jsq,1)=sigma1;
                    param(jsq,2)=t;
                    param(jsq,3)=h1;
                    profile=profile1;
                    bgdata=bgdata1;
                    bgvalue(jsq,:)=bgdata1;
                    profiledata(jsq,:)=profile1;
                else
                    param(jsq,1)=sigma;
                    param(jsq,2)=t;
                    param(jsq,3)=h;
                    bgvalue(jsq,:)=bgdata;
                    profiledata(jsq,:)=profile;
                end
            else
                    sigma=sigma1;
                    h=h1;
                    param(jsq,1)=sigma1;
                    param(jsq,2)=t;
                    param(jsq,3)=h1;
                    profile=profile1;
                    bgdata=bgdata1;
                    bgvalue(jsq,:)=bgdata1;
                    profiledata(jsq,:)=profile1;
            end

        else num==1 %radius 9nm/pixel
            jsq=jsq+1;
            t1=((th-tl)*C(jsq+500)+tl)/13;
            mu_est=50*sigma;
            gauss_psf=norm_density(x, mu_est, sigma);
            exp_decay= zeros(length(x),1); 
            exp_decay(1:round(50*sigma))=0;
            exp_decay(round(50*sigma):end)=exp(-(x(round(50*sigma):end)-round(50*sigma))/t1);
            exp_gauss=conv(exp_decay,gauss_psf);
            [h1,profile1,bgdata1]=sqsum(exp_gauss,profilemeanI,l1,l2,center_value,bgvalue_record,maxI);
            Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
            if Z<1
                num1=D(jsq+500);
                if Z>num1
                    t=t1;
                    h=h1;
                    param(jsq,1)=sigma;
                    param(jsq,2)=t1;
                    param(jsq,3)=h1;
                    profile=profile1;
                    bgdata=bgdata1;
                    bgvalue(jsq,:)=bgdata1;
                    profiledata(jsq,:)=profile1;
                else
                    param(jsq,1)=sigma;
                    param(jsq,2)=t;
                    param(jsq,3)=h;
                    bgvalue(jsq,:)=bgdata;
                    profiledata(jsq,:)=profile;
                end
            else
                    t=t1;
                    h=h1;
                    param(jsq,1)=sigma;
                    param(jsq,2)=t1;
                    param(jsq,3)=h1;
                    profile=profile1;
                    bgdata=bgdata1;
                    bgvalue(jsq,:)=bgdata1;
                    profiledata(jsq,:)=profile1;
            end

        end
      end
        
end

