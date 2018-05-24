%%%% 2nd program
syn_syn=fitsread('syn_syn_1024_with_monopole.fits','table');
dust_dust=fitsread('dust_dust_1024_with_monopole.fits','table');
%dust_syn=fitsread('dust_syn_haslam_with_monopole.fits','table');
%syn_syn_TT=fitsread('syn_syn_TT_with_monopole.fits','table');
%dust_dust_TT=fitsread('dust_dust_TT_with_monopole.fits','table');
%dust_syn_TT=fitsread('dust_syn_TT_with_monopole.fits','table');
dust_dust{1}=dust_dust{1};
dust_syn=sqrt(syn_syn{1}.*dust_dust{1});


%%%%%%%%% planck syn_pol_1024 is needed!!


F(2,2)=0;
F(1,1)=syn_syn{1}(1)/4/pi;
F(2,2)=dust_dust{1}(1)/4/pi;
F(1,2)=dust_syn(1)/4/pi; 
F(2,1)=F(1,2);

scale=[(70/30)^beta_s 0; 0 (70/353)^(beta_d+1)*(exp(353/v_Td)-1)/(exp(70/v_Td)-1)];
F_scale=scale*F*scale;
F_scale=[0.0025,0.082;0.082,3.20]
%F=scale^-1*F_scale*scale^-1
ell=2:512;


%figure(fg1)
hold on
loglog(ell,dust_dust{1}(3:end)'.*ell.*(ell+1)/2/pi*((65/353)^(beta_d+1)*(exp(353/v_Td)-1)/(exp(65/v_Td)-1))^2,'--')

loglog(ell,syn_syn{1}(3:end)'.*ell.*(ell+1)/2/pi*(65/30)^(2*beta_s),':')

loglog(ell,dust_syn(3:end)'.*ell.*(ell+1)/2/pi*(65/30)^(beta_s)*((65/353)^(beta_d+1)*(exp(353/v_Td)-1)/(exp(65/v_Td)-1)),'-.')
