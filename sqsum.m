function [h_m,profile_m,background_m]=sqsum(exp_gauss,profilemeanI,l1,l2,centergaus,bgvalue_record,maxI)
   [pImax,pIexp]=max(exp_gauss);
   lmin=min(l1,l2);% 13nm/pixel
   lmax=max(l1,l2);%13nm/pixel
   mask=ones(5,1);%13nm/pixel to 65nm/pixel
%    [pmax,pI]=max(profilemeanI);
   pI=round(centergaus);
   rng('shuffle');
   b_frac=rand(4000,10);
   profilemeanI_n=profilemeanI/max(profilemeanI);
   exp_data=exp_gauss(pIexp-round(lmin*5.5):pIexp+round(lmax*2.6));
   if (pI+round(2.3*lmax/5))<length(profilemeanI) && (pI-round(3.2*lmin/5))>1
   profile_data= profilemeanI_n(pI-round(3.2*lmin/5):pI+round(2.3*lmax/5));
   elseif (pI-round(3.2*lmin/5))>0
   profile_data= profilemeanI_n(pI-round(3.2*lmin/5):end);
   end
   lw=length(exp_data)-length(profile_data)*5;
    
    for j=0:1:lw
        IAx=((j+1):1:(j+length(profile_data)*5));
        expdata=exp_data(IAx);
        exp_gauss_c=conv(expdata,mask,'valid');
        eg_c=exp_gauss_c(1:5:end);
        egc=eg_c/max(eg_c);
        %% add poisson noise to data
        eg_data=egc*maxI;
        alpha=0.75; %photon number
            for ip=1:1:40
            %% add photon shot noise
            avg_shot=poissrnd(eg_data,1,length(eg_data));
        %% alpha conversion factor
            avg_read=uint16((avg_shot*alpha)/0.53);% change to ADU unit
            [emax,mind]=max(avg_read);
            egcn=double(avg_read)/double(emax);
            %% choose two bg for two tails
            bg_im1=0.15*(b_frac(1001+j:3000+j,3)-0.5)+bgvalue_record;
            bg_im2=0.15*(b_frac(1001+j:3000+j,6)-0.5)+bgvalue_record;
            dnum=length(bg_im1);
            eg_rep=repmat(egcn,[dnum,1]);
            bg_im_mat1=repmat(bg_im1,[1,mind]);
            bg_im_mat2=repmat(bg_im2,[1,length(egcn)-mind+1]);
            imag1=eg_rep(:,1:mind)+bg_im_mat1;
            imag2=eg_rep(:,mind:length(egcn))+bg_im_mat2;
            [maxline1,line_layer1]=max(imag1');
            [maxline2,line_layer2]=max(imag2');
            imagenorm1=imag1./maxline1';
            imagenorm2=imag2./maxline2';
            bg1=bg_im_mat1./maxline1';
            bg2=bg_im_mat2./maxline2';

            expline_norm1=repmat(profile_data(1:mind),[dnum,1]);
            expline_norm2=repmat(profile_data(mind:length(egc)),[dnum,1]);

            diff_image1=(imagenorm1-expline_norm1).^2;
            diff_image2=(imagenorm2-expline_norm2).^2;

            diff_mid1=(imagenorm1-expline_norm1).^2.*abs(expline_norm1)*5;
            diff_mid2=(imagenorm2-expline_norm2).^2.*abs(expline_norm2)*5;
%             diff_image=[diff_image1,diff_image2];

            diff_image_s1=sum(diff_mid1')+sum(diff_image1');
            diff_image_s2=sum(diff_mid2')+sum(diff_image2');

            [minvalue1,minind1]=min(diff_image_s1);
            [minvalue2,minind2]=min(diff_image_s2);

            delta(ip)=min(diff_image_s1)+min(diff_image_s2);
            minbg(ip,1:mind)=imagenorm1(minind1,:);
            minbg(ip,mind:length(egc))=imagenorm2(minind2,:);

            bgvalue(ip,1)=bg1(minind1(1),1);
            bgvalue(ip,2)=bg2(minind2(1),1);
            end
           
     [h(j+1),expind]=min(delta);
     profile(j+1,:)=minbg(expind,:);
     background(j+1,:)=bgvalue(expind,:);
    
    end
    
    [h_m,expind_m]=min(h);
    profile_m=profile(expind_m,:);
    background_m=background(expind_m,:);
    
end