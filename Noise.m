% 1st Program
clear;
fre=[45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 375, 435, 555, 675, 795]; %in GHz
Del_fre=[15 15 15 15 15 15 15 15 15 15 15 15 195 195 195]; %in GHz
n_det=[64 300 400 550 750 1150 1800 575 375 100 64 64 64 64 64];
FWHM=[23.3, 14.0, 10.0, 7.8, 6.4, 5.4, 4.7, 4.1, 3.7, 3.3, 2.8, 2.4, 1.9, 1.6, 1.3]; %in arcmin
FWHM_rad=FWHM/60/180*pi; 
RJ_temp=[4.98 2.36 2.03 1.68 1.38 1.07 0.82 1.40 1.70 3.25 4.05 4.12 1.23 1.28 1.31]; %in muK*arcmin
RJ_pol=[8.61, 4.09, 3.50, 2.90, 2.38, 1.84, 1.42, 2.43, 2.94, 5.62, 7.01, 7.12, 3.39, 3.52, 3.60];%in muK*arcmin
RJ_temp_rad=RJ_temp/60/180*pi;  %%now in muK*rad
RJ_pol_rad=RJ_pol/60/180*pi;  %%now in muK*rad
sigma_temp_weighted=[5.25 2.73 2.68 2.63 2.67 2.63 2.64 6.08 10.1 26.9 68.6 149 227 1320 8070];
sigma_pol_weighted=[9.07 4.72 4.63 4.55 4.61 4.54 4.57 10.5 17.4 46.6 119 258 626 3640 22200];
sigma_pol_weighted_rad=sigma_pol_weighted/60/180*pi;    %%now in muK*rad
sdet=30; %in muK*s^1/2
nside=512;
fsky=0.7;



lmin=2;
lmax=512;
ell=lmin:lmax;



%[23.3, 14.0, 10.0, 7.8, 6.4, 5.4, 4.7, 4.1, 3.7, 3.3, 2.8, 2.4, 1.9, 1.6, 1.3]
%[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 195, 195, 195]
%[45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 375, 435, 555, 675, 795]
%v=0:0.01:1000;

v_ref_syn_pol=30; % Planck 2015 reference frequency of synchrotron polarization in GHz
v_ref_dust_pol=353; % Planck 2015 reference frequency of dust polarization in GHz
v_ref=150;
beta_s=-3.1;
%beta_s=-3; %Errard2011 used -3
beta_d=1.59;
%beta_d=1.65; %Errard2011 used -1.65
v_Td=407.9638; % 19.6K
%v_Td=407.9638*18/19.6;  %Errard2011 used 18K
v_cmb=56.729; % 2.72548K

%A_cmb = @(v) (v/v_ref).^2.*exp((v-v_ref)/v_cmb)./(exp(v/v_cmb)-1).^2*(exp(v_ref/v_cmb)-1)^2;
A_cmb = @(v) (v/v_cmb).^2.*exp(v/v_cmb)./(exp(v/v_cmb)-1).^2;
A_sync = @(v) (v/v_ref_syn_pol).^beta_s;
A_dust = @(v) (v/v_ref_dust_pol).^(beta_d+1)*(exp(v_ref_dust_pol/v_Td)-1)./(exp(v/v_Td)-1);


A(length(fre),3)=0; % 1:cmb, 2:sync, 3: dust
for i=1:length(fre)
        A(i,1)= integral(A_cmb,fre(i)-Del_fre(i)/2,fre(i)+Del_fre(i)/2)/Del_fre(i);
        A(i,2)= integral(A_sync,fre(i)-Del_fre(i)/2,fre(i)+Del_fre(i)/2)/Del_fre(i);
        A(i,3)= integral(A_dust,fre(i)-Del_fre(i)/2,fre(i)+Del_fre(i)/2)/Del_fre(i);
end

%%%Noise post separation
N_noise(length(fre),length(fre),length(ell))=0;
for j=1:length(ell)
  for i=1:length(fre)
    N_noise(i,i,j)=RJ_pol_rad(i)^2*exp((j+lmin)*(j+lmin)*FWHM_rad(i)^2/8/log(2));
  end
end

N_post_l(1:length(ell))=0;
for i=1:length(ell)
  something=(A'*N_noise(:,:,i)^-1*A)^-1;
  N_post_l(i)=something(1,1);
end
%N_pre=1/sum(RJ_pol_rad.^-2);

N(length(fre),length(fre))=0;
for i=1:length(fre)
%N(i,i)=4*(30)^2*12*nside^2/4/365/24/3600/n_det(i);
N(i,i)=(RJ_pol(i))^2/41252.96*12*nside^2/60/60;
end

N_degrade=(A'*N^-1*A)^-1;
N_inv_weighted=1/(A(:,1)'*N^-1*A(:,1)); %%% same as N_pre
N_pre=1/sum(sigma_pol_weighted.^-2*41252.96/12/nside^2*60*60);  %%%in muK
degradation=N_degrade(1,1)/N_pre    

N_vv(length(fre),length(ell))=0;
for j=1:length(ell)
  for i=1:length(fre)
    N_vv(i,j)=sigma_pol_weighted_rad(i)^(2)*exp((j+lmin)*(j+lmin)*FWHM_rad(i)^2/8/log(2));
  end
end


%%%End Noise post separation

%%%%
N_ell(1:length(ell))=0;
for j=1:length(ell)
   N_ell(j)=1/sum(N_vv(:,j).^-1);
end



N_ell=degradation*N_ell;
Nellell=N_ell.*ell.*(ell+1)/2/pi;
%fg1=figure;
loglog(ell,Nellell,'r')
hold on
%loglog(ell,N_post_l.*ell.*(ell+1)/2/pi)
%%%% So N_post_l and N_ell here are the same