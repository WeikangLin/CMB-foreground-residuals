%%%%%% 4th program, caculate the foreground residuals

alpha(3,3,2)=0;
alpha(:,:,1)=-N_degrade*A'*N^-1*A_beta(:,:,1);
alpha(:,:,2)=-N_degrade*A'*N^-1*A_beta(:,:,2);


%%% The 2*2*2*2 matrix Kapper
Kapper(2,2,2,2)=0;
for i=1:2
    for j=1:2
        for m=1:2
            for n=1:2
                Kapper(i,j,m,n)=alpha(1,m+1,i)*alpha(1,n+1,j);
            end
        end
    end
end

%%%%%% calculating the foreground residuals
C_fg_res(length(ell))=0;
C_dust_res(length(ell))=0;
C_syn_res((length(ell)))=0;
C_dust_syn_res((length(ell)))=0;
for i=1:length(ell)
    C_syn_res(i)=Sigma(1,1)*(Kapper(1,1,1,1)*syn_syn{1}(i+2)+2*Kapper(1,1,1,2)*dust_syn(i+2)+Kapper(1,1,2,2)*dust_dust{1}(i+2));
    C_fg_res(i)=C_syn_res(i);
    C_dust_syn_res(i)=2*Sigma(1,2)*(Kapper(1,2,1,1)*syn_syn{1}(i+2)+2*Kapper(1,2,1,2)*dust_syn(i+2)+Kapper(1,2,2,2)*dust_dust{1}(i+2));
    C_fg_res(i)=C_fg_res(i)+C_dust_syn_res(i);
    C_dust_res(i)=Sigma(2,2)*(Kapper(2,2,1,1)*syn_syn{1}(i+2)+2*Kapper(2,2,1,2)*dust_syn(i+2)+Kapper(2,2,2,2)*dust_dust{1}(i+2));
    C_fg_res(i)=C_fg_res(i)+C_dust_res(i);
end

%loglog(ell,C_dust_res.*ell.*(ell+1)/2/pi,'r')
%loglog(ell,C_dust_syn_res.*ell.*(ell+1)/2/pi)
%loglog(ell,C_syn_res.*ell.*(ell+1)/2/pi,'g')

loglog(ell,C_fg_res.*ell.*(ell+1)/2/pi,'b')
hold on
total=N_post_l'+C_fg_res';
loglog(ell,total'.*ell.*(ell+1)/2/pi,'k','LineWidth',1.25)

save('total_BB.txt','total','-ascii')

xlabel('l','FontSize', 18)
ylabel('l(l+1)C_l/2\pi (\muK^2)','FontSize', 18)
axis([2,500,1E-7,1])
grid on

axis_x_y=axis;
text(0.3*axis_x_y(2),1E-6,'COrE','FontSize', 18)
set(gca,'FontSize',12)
ax=gca;
ax.XTick = [2 10 100 500];
%ax.XRuler.MinorTick = [3 4 5 6 7 8 9 20 30 40 50 60 70 80 90 200 300 400];
%ax.YRuler.MinorTick = [1e-5 1e-3 1e-1];
