%% ========================== FIGURES =====================================
% Output for the presentatoins and paper
% =========================================================================
% modified on October 10 2018
%
%% ************************** homeostatic optimization results ************
%---------------------------- mass fraction--------------------------------
h1=figure;
hold on;
for i=1:Newgen
    [e, c, m] = mass_fracs(2*Radius(i));
    semilogx(2*Radius(i),e,'ob');
    semilogx(2*Radius(i),m,'or');
    semilogx(2*Radius(i),1 - e - m,'om');
    grid on;
end
xlabel('log(D[\mum])'); ylabel('Mass fraction');
legend('elastin','SMC','collagen');
hold off;
% FigModify(h1,'MassFraction')

% for generations 
k = 1:1:Newgen;  
%------------------------h/2R ratio --------------------------------------
h2=figure;
plot(2*Radius(k)*100,ratio(k),'o-','LineWidth',2);
set(gca,'XDir','reverse');xlabel({'2R (cm)'},'FontSize',12);
ylabel({'h/2R'},'FontSize',12);hold on;
legend({'Thickness-to-diameter ratio'},'FontSize',12);
grid on;
% FigModify(h2,'ratioHto2R')

%------------------------Lenght vs gen-------------------------------------
h3=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k, Length(k).*100,'o-','LineWidth',2);   % factor 0.5 from Olufsen
plot(k, Length(k).*2*100,'--','LineWidth',2); % back to Olufsen curve
plot(1, L0*100,'*','LineWidth',2); 
plot(1, 1.807,'*','LineWidth',2); 
legend({'Approximation 0.5L_{Oluf}','L_{Oluf}','RIA [Olufsen-2012]','12th order [Huang-1993])'},'FontSize',12);
xlabel({'gen.'},'FontSize',12);ylabel({'Length (cm)'},'FontSize',12);
grid on;
% FigModify(h3,'LvsGen')

%----------------------- area ratio----------------------------------------
%(A1+A2)/A0=2A1/A0=2R1^2/R0^2
h4=figure; 
set(gca,'xscale'); xlim([0 Newgen]); hold on;
eta01=1.2; eta02=1.3;
eta=ones(1,Newgen);
for m=2:Newgen
    eta(m)=2*Radius(m)^2/Radius(m-1)^2;
end
etaR1(1:Newgen)=eta01;etaR2(1:Newgen)=eta02;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(2:Newgen, eta(2:Newgen),'o--','LineWidth',2); 
plot(2:Newgen,etaR1(2:Newgen),'r--','LineWidth',2);
plot(2:Newgen,etaR2(2:Newgen),'r--','LineWidth',2);
legend({'Optimization Results','open-end: 1.2< area ratio <1.3 [Hollander-2001]'},...
    'FontSize',12);
xlabel({'gen.'},'FontSize',12);ylabel({'Daughter-to-parent area ratio'},'FontSize',12);
axis([1 Newgen 1.1 1.4]);
grid on;
% FigModify(h4,'AreaRatio')

%-----------------------Eh/R0----------------------------------------------
h5=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
EhrKall(1:Newgen)=EhrK; EhrYall(1:Newgen)=EhrY; EhrQall(1:Newgen)=EhrQ;
plot(k, Thickness(k).*YoungMod_tt(k)./RzeroP(k)./1000,'o--','LineWidth',2); 
plot(k, EhrKall./1000,'--','LineWidth',2); 
plot(k, EhrYall./1000,'-.','LineWidth',2);
errorbar(1,EhrEx_mean/1000,abs(min(EhrEx_err/1000)),abs(max(EhrEx_err/1000)),...
    'o','LineWidth',2);
legend({'Optimization Results','Eh/R0 [Krenz-2003]','Eh/R0 [Yen-1990]',...
    'E_{\theta\theta}h/R0 (MSU Experiments)'},'Location','NorthEast',...
    'FontSize',12);
xlabel('gen.','FontSize',12);ylabel('E_{\theta\theta}h/R0 (kPa)',...
    'FontSize',12);
grid on;
axis([0 Newgen 2 14]);

%---------------homeostatic values vs radius (and generations)-------------
h6=figure;
subplot(3,2,1);
%     set(gca,'xscale');
set(gca,'XDir','reverse');
xlabel('d (cm)');ylabel('E_{\theta\theta} (kPa)');hold on;
plot(2*Radius(k)*100,YoungMod_tt(k)/1000,'o-'); %,'LineWidth',1.5);

subplot(3,2,2); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('\sigma_h (kPa)');hold on;
plot(2*Radius(k)*100,sigma_h(k)/1000,'o-'); %,'LineWidth',1.5);

subplot(3,2,3); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('\tau (Pa)');hold on;
plot(2*Radius(k)*100,shear(k),'o-'); %,'LineWidth',1.5);

subplot(3,2,4); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('Pmid (mmHg)'); hold on;
plot(2*Radius(k)*100,Pmid(k)/133.32,'o-'); %,'LineWidth',1.5);

subplot(3,2,5); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('h/d');hold on;
plot(2*Radius(k)*100,ratio(k),'o-'); %,'LineWidth',1.5);

subplot(3,2,6); set(gca,'xscale');xlabel('gen.');ylabel(' Radius exponent \xi');hold on;
plot(1:Newgen-1,ksi(1:Newgen-1),'o-'); %,'LineWidth',1.5);
xlim([0 Newgen]);

% suptitle(['No. generations = ',num2str(Newgen),' for D_{min} = ',num2str(2*Rmin*100),'cm']);
hold off;

%----------------------geometrical values vs generations-------------------
h7=figure;

subplot(2,1,1); set(gca,'xscale');  xlim([0 Newgen]); hold on;
plot(k, 2*Radius(k).*100,'o-'); 
plot(1 ,2*R0*100,'*');
plot(1 ,0.271,'*');
legend('homeostatic values','2R_0 (RIA, Olufsen)','2R (12 order, Huang)');
xlabel('gen.');ylabel('d (cm)');

subplot(2,1,2); set(gca,'xscale');  xlim([0 Newgen]); hold on;
plot(k, Thickness(k).*100,'o-'); 
plot(1 ,H0*100,'*');
legend('homeostatic values','H_0=0.1*R_0 (Banks)');
xlabel('gen.');ylabel('h (cm)');

% suptitle(['No. generations = ',num2str(Newgen),' for D_{min} = ',num2str(2*Rmin*100),'cm']);
hold off;

%--------------------material properties vs generations--------------------
h8=figure;

subplot(3,1,1); set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k, YoungMod_tt(k)./1000,'o-'); 
plot(k, YoungMod_zz(k)./1000,'o-'); 
plot(1, E0/1000,'*'); 
legend('homeostatic values-E_{\theta\theta}','homeostatic values-E_{zz}','? E_{\theta\theta} - est.(3/4/5)*C_{\theta\theta}_{AA-Roccabianca}');
xlabel('gen.');ylabel('Young modulus (kPa)');

subplot(3,1,2); set(gca,'xscale'); xlim([0 Newgen]); hold on;
nu0(1:Newgen)=0.5;
plot(k, nu_tz(k),'o-'); 
plot(k, nu_zt(k),'o-'); 
plot(k, nu0,'--'); 
legend('nu_{\theta z}','nu_{z \theta}','nu = 0.5 isotropic incompressible');
xlabel('gen.');ylabel('Poisson ratio ');

subplot(3,1,3); set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k, Thickness(k).*YoungMod_tt(k)./1000,'o--'); 
plot(1, H0*E0/1000,'*'); 
legend('homeostatic values','? E_{\theta\theta}(est.)*H0','Location','NorthEast');
xlabel('gen.');ylabel('E_{\theta\theta}*h (kPa*m)');

hold off;

%% ************************** pulsatile hemodynamics results ************
%------------------------------PWV-----------------------------------------
h9=figure;
set(gca,'xscale');  xlim([0 Newgen]); hold on;
plot(k, c_R2,'o-','LineWidth',2); 
plot(k, c_R10,'o-','LineWidth',2); 
plot(k ,c0_MK,'.-k','LineWidth',2);
% plot(k ,c0_tt,'--k','LineWidth',2);
plot(1 ,PwvB,'*','LineWidth',2);
plot(1 ,PwvM,'*','LineWidth',2);
legend({'at 1st mode \omega=2\pi/T','at 9th mode \omega=2\pi9/T',...
 'c^{MK} with E_{\theta\theta}','[Banks-1978]',...
 '[Milnor-1969]'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('Pulse wave velocity (m/s)','FontSize',12);
grid on;
axis([0 Newgen 0 3]);

%------------------------------Delta---------------------------------------
h10=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k,delta1(k),'o-','LineWidth',2);
plot(k,delta9(k),'o-','LineWidth',2);
% legend({'delta=max(longVel)/c_{r1}','deltaMK=max(longVel)/c_{MK}'},'Location','NorthEast','FontSize',12);
legend({'{\delta} at 1st mode','{\delta} at 9th mode'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('{\delta}=max(v^f_{z})/c','FontSize',12);
grid on;

%--------------------------Womersley number--------------------------------
h11=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k,Alpha(2,k),'o-','LineWidth',2);
plot(k,Alpha(10,k),'o-','LineWidth',2);
legend({'{\alpha} at 1st mode  \omega=2\pi/T','{\alpha} at 9th mode \omega=2\pi9/T'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('Womersley number','FontSize',12);
grid on;


% time points
t = linspace(0,T-T/Nt,Nt);
%--------------------------Flow and Pressure-------------------------------
h12=figure;
col=jet(Newgen); 

subplot(2,1,1);
hold on;
for k=1:Newgen
    if k==1 
        plot(t,qInpTime(:,k)*10^6,'color',col(k,:),'LineWidth',2);
    else
        plot(t,qInpTime(:,k)*10^6,'color',col(k,:));
    end
    hold on
end
title(['Input Flow over ',num2str(Newgen),' gen.'],'FontSize',12);
xlabel('Time (s)','FontSize',12); ylabel('Flow (cm^3/s)','FontSize',12); grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T -5 50]);

subplot(2,1,2);
hold on;
for k=1:Newgen
    if k==1
        plot(t,pTermTime(:,k)/133.32,'color',col(k,:),'LineWidth',2);    
    else
        plot(t,pTermTime(:,k)/133.32,'color',col(k,:));
    end
    hold on
end
hold off
title(['Terminal Pressure over ',num2str(Newgen),' gen.'],'FontSize',12);
xlabel('Time (s)','FontSize',12);  ylabel('Presure (mmHg)','FontSize',12);   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T 5 25]);

%------------------------ Root pressure------------------------------------
h13=figure;
plot(t,pInpTime(:,1)/133.32,'-','LineWidth',2);
hold on
plot(t,pTermTime(:,1)/133.32,'-','LineWidth',2);
legend({'Input p at 1st gen.','Terminal p at 1st gen.'},'FontSize',12);
xlabel('Time (s)','FontSize',12);  ylabel('Presure (mmHg)','FontSize',12);   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T 5 25]);

%------------------------ Root flow----------------------------------------
h14=figure;
plot(t,qInpTime(:,1)*10^6,'-','LineWidth',2);
hold on
plot(t,qTermTime(:,1)*10^6,'-','LineWidth',2);
legend({'Input q at 1st gen.','Terminal q at 1st gen.'},'FontSize',12);
xlabel('Time (s)','FontSize',12);  ylabel('Flow (cm^3/s)','FontSize',12);   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T -5 50]);

%---------------------- Impedance in time domain---------------------------
h15=figure;
% plot(t,z./(133.32*10^6),'--','LineWidth',2);  hold on
plot(t,zinp_total./(133.32*10^6),'-','LineWidth',2);  hold on
plot(t,zterm_total./(133.32*10^6),'-.','LineWidth',1);  hold on
plot(t,zchar(:,1)./(133.32*10^6),'--','LineWidth',1);  hold off
grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
legend({'Z_{inp}','Z_{term}','Z_c-(pulsatile)'},'FontSize',12);
xlabel('Time (s)','FontSize',12); ylabel('Impedance(mmHg*s/cm^3)','FontSize',12);
axis([0 T -4 10]);

%---------------------- Impedance in frequency domain----------------------
h16=figure;
kk = 0:1:NumModes-1;
% suptitle('Input impedance in frequency domain');
subplot(1,2,2); plot(kk,angle(Zn),'o','LineWidth',2); grid on;
% set(gca,'ytick',-pi:pi/4:pi); 
% set(gca,'yticklabel',{'-\pi','-\pi/2','-\pi/4','0','\pi/4','\pi/2','\pi'}); 
xlabel('Frequency modes','FontSize',12); ylabel('Phase Angle','FontSize',12);
subplot(1,2,1); plot(kk,abs(Zn)./(133.32*10^6),'o','LineWidth',2); grid on;
set(gca,'xtick'); xlabel('Frequency modes','FontSize',12); ylabel('Modulus (mmHg*s/cm^3)','FontSize',12);
