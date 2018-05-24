%%%%%% 3rd program
%%%% using N 
per=0.0001;


%%%%%%This is extra to check N
%N(length(fre),length(fre))=0;
%for i=1:length(fre)
%N(i,i)=4*(30)^2*12*nside^2/4/365/24/3600/n_det(i);
%end
%N_degrade=(A'*N^-1*A)^-1;
%%%%%%This is extra to check N


A_beta(length(fre),3,2)=0;  %derivative of A 1 for syn, 2 for dust

%for syn
A_beta(:,:,1)=A;
A_sync_p = @(v) (v/v_ref_syn_pol).^(beta_s*(1+per));
for i=1:length(fre)
   A_beta(i,2,1)= integral(A_sync_p,fre(i)-Del_fre(i)/2,fre(i)+Del_fre(i)/2)/Del_fre(i);
end
A_beta(:,:,1)=(A_beta(:,:,1)-A)/per/beta_s;

%for dust
A_beta(:,:,2)=A;
A_dust_p = @(v) (v/v_ref_dust_pol).^(beta_d*(1+per)+1)*(exp(v_ref_dust_pol/v_Td)-1)./(exp(v/v_Td)-1);
for i=1:length(fre)
   A_beta(i,3,2)= integral(A_dust_p,fre(i)-Del_fre(i)/2,fre(i)+Del_fre(i)/2)/Del_fre(i);
end
A_beta(:,:,2)=(A_beta(:,:,2)-A)/per/beta_d;

%%%constructing the long matrix, only 2nd and 3rd components are non-zero
Massive(3,3,2,2)=0;
Massive(:,:,1,1)=A_beta(:,:,1)'*N^-1*A*N_degrade*A'*N^-1*A_beta(:,:,1)-A_beta(:,:,1)'*N^-1*A_beta(:,:,1);
Massive(:,:,1,2)=A_beta(:,:,1)'*N^-1*A*N_degrade*A'*N^-1*A_beta(:,:,2)-A_beta(:,:,1)'*N^-1*A_beta(:,:,2);
Massive(:,:,2,1)=Massive(:,:,1,2)';
Massive(:,:,2,2)=A_beta(:,:,2)'*N^-1*A*N_degrade*A'*N^-1*A_beta(:,:,2)-A_beta(:,:,2)'*N^-1*A_beta(:,:,2);

Reduce=Massive(2:3,2:3,:,:);%%%%%%% since the 1st component of Massive matrix is 0's

%%%% F matrix%%%%%%%%%%%%%%%%%%% Need to use a new one for our calculation
%F_P06_70=[0.0025 0.082; 0.082 3.2];
%F_me=[0.1495 0.3729; 0.3729 2.0683];
F_dif_ref=F;
%%rescaling

scale70=[(v_ref_syn_pol/70)^beta_s 0; 0 (v_ref_dust_pol/70)^(beta_d+1)*(exp(70/v_Td)-1)/(exp(v_ref_dust_pol/v_Td)-1)];
scale150=[(v_ref_syn_pol/150)^beta_s 0; 0 (v_ref_dust_pol/150)^(beta_d+1)*(exp(150/v_Td)-1)/(exp(v_ref_dust_pol/v_Td)-1)];
F_scale150=scale150*F*scale150';


%%%%Sigma_inv and Sigma matrix
Sigma_inv(2,2)=0;
Sigma_inv(1,1)=-1*trace(Reduce(:,:,1,1)*F_dif_ref);
Sigma_inv(1,2)=-1*trace(Reduce(:,:,1,2)*F_dif_ref);
Sigma_inv(2,1)=Sigma_inv(1,2);
Sigma_inv(2,2)=-1*trace(Reduce(:,:,2,2)*F_dif_ref);

Sigma_inv=Sigma_inv*fsky*12*nside^2; %%%%Since I am following Errard2011

Sigma=Sigma_inv^(-1);

d_beta_s=Sigma(1,1)^0.5
d_beta_d=Sigma(2,2)^0.5

x=d_beta_s;
y=d_beta_d;
xy=Sigma(1,2);

a=sqrt((x^2+y^2)/2+sqrt((x^2-y^2)^2/4+xy^2));
b=sqrt((x^2+y^2)/2-sqrt((x^2-y^2)^2/4+xy^2));
theta=atan(2*xy/(x^2-y^2))/2;


%t=linspace(0,2*pi);
%fh2 = figure;
%plot(beta_s+1.52*a*cos(t)*cos(theta)-1.52*b*sin(t)*sin(theta),beta_d+1.52*b*sin(t)*cos(theta)+1.52*a*cos(t)*sin(theta))
%hold on
%plot(beta_s+2.48*a*cos(t)*cos(theta)-2.48*b*sin(t)*sin(theta),beta_d+2.48*b*sin(t)*cos(theta)+2.48*a*cos(t)*sin(theta))
%plot(nu_0+3.44*a*cos(t)*cos(theta)-3.44*b*sin(t)*sin(theta),r+3.44*b*sin(t)*cos(theta)+3.44*a*cos(t)*sin(theta))


