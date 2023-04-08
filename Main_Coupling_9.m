clear; close all; clc;
Const = load('Const_Mass.mat');
% Homeostatic Optimization and Hemodynamics in arterial tree
%   - Extension of Murray's Law for homeostasis (steady flow)
%   - Symmetric bifurcating tree
%   - Pulsitile flow postprocessing as for Deformable wall Womersley theory 
%       (longitudinally constrained)
%   - Axisymmetric straight cylindrical vessels

% code by Hamid and Vasilina on August 2 2018

% Important Note: 
%   - System of Units (m,s,kg)
%   - mmHg = 1,333.2 dyn/cm2 = 1,33.32 Pa

%% Wall Tissue Parameters
% Assumptions: 
%   1) we have 4 collagen fiber families in 0, 90, 45, and -45 degrees. 
%   2) SMCs that are circumferentially oriented. 
%   3) isotropic elastin matrix.
%   4) The mass fraction in generations is known.
%   5) Ratio of collagen to SMC is constant throughout the arterial tree. 
%   6) The deposition stretch of SMC and Collagen fibers are constant.
%   8) Stiffness components are computed using Small on Large theory (SoL).
%   9) The metabolic cost of maintenance of collagen and SMCs is the same 
%   and constant. 
%   7) Blood viscosity and metabolic cost (of ?) is constant. 

% set some global parameters (NOT variable!):
global H0 rho_w beta gamma R0 Rmin alpha_t mu lmax lmin S ...
       Ghe1 Ghe2 Ghc Ghm fiber_ang

% Metabolic costs
alpha_t = 2*746;    % W/m^3, Taber 1998 (Biophys J.), porcine carotid artery
beta = 160;         % W/m^3, Lindstrom et al (2014)
gamma = 0.00891267681;    % J*s/m^3, Lindstrom et al (2014)

% Active wall
S =  2.0e+004;
lmax = 1.2;     %  nondim., maximumu tension in active tone
lmin = 0.7;     %  nondim., minimum tension in active tone

% Wall density
rho_w = 1060;   % kg/m3 

% Mechanical_properties (see function mechanical_properties_PA)
Ghe1 = 1.27;
Ghe2 = 1.27;
fiber_ang = 45*pi/180; % angle of helical collagen fibers
Ghc = 1.154; %HG changed from 1.034. June 22nd
Ghm = 1.21;

%% Tree Parameters
% Reduction in elastin's mass fraction in each generation
f = 0.5; 

% initial radius of root vessel [Olufsen 2012]RIA: 0.011/2; 
R0 = 0.0055; %m
% assumed initial thickness of root vessel (as disscussed with Baek)
H0 = 0.07*R0;

% Minimum vesel radius to define terminal vessel (or number of generations)
% [Olufsen 2012]: Rmin=50microns=0.005cm=0.00005mn 
%   - arterioles, not cappilaries 
Rmin = 0.00005; % m

% preliminary maximum number of generations
% >=19
N_gen = 19;

%% Reference Values for Results Comparisons
L0 = 0.0125; %[Olufsen 2012]RIA

% Stiffness parameter E11*H*Rp0 = 3/4/Lambda for isotropic incompressible
%   - Lanbda - distensibility parameter lambda
%   - Rp0 - radius at zero pressure
% (Krenz 2003) Lambda=0.02/mmHg
% (Yen1990 according to Qureshi 2014) lambda=0.012 mmHg
%   - HG: average lambda from Yen et al (1990) 0.0146 mmHg
EhrY = 3/(4*0.012)*133.32;

% Ehr from porcine experiments
EhrExperiment = [10530.92136,9496.887786,9311.097003,8679.807346,8432.782151]; 
EhrEx_mean = mean(EhrExperiment);EhrEx_err = EhrExperiment - EhrEx_mean;

% questionable assumptions
EhrK = 37.5*133.32; %Olufsen 2012, Eh/r0=3/4/Lambda, 
EhrQ = 195*133.32;  %Qureshi 2014 

% Elastic modulus at pulmonary artery
%   as 3C/4/5=0.1125MPa=112.5kPa (from C=0.75MPa from Roccabiance ATA)
%   assuming poisson's ratio=0.5
E0 = 112500;

% wave speed - pulse wave velocity
% Banks (1978): Moens-Korteweg PVW for 20-30yo humans c0=2.24m/s
% Milnor (1969): c0=1.68m/s;
PwvB = 2.24; PwvM = 1.68;

% Normal mean arterial pressure 15, diastolic 10, systolic 25 mmHg;
Pm0=15*133.32;  Pd0=10*133.32; Ps0=25*133.32;

%% Hemodynamics Parameters
% Assumptions for Pulsatile hemodynamics:
%   - Find pulsatile solution from Womersley theory (longitudinally 
%     constrained wall)
%   - Given: input flowaveform and zero oscillatory terminal pressure at
%     the last generation (open-end);
%   - Steps:
%   1) Knowing geometry (R,L,symmetry)and wall properties (E,h or Stiffness 
%   matrix component) compute impedances in each segment;
%   2) Knowing impedance and steady solution, compute the total p and q for
%   each segment;
%   3) Postprocess results
%   4) Do that for the realistic pulmonary vessels parameters;
%   - Terminal reflection coefficient GammaVal
%       0   - matching impedance (Z0=ZT=Zinp), 
%       1   - closed end (ZT>>Z0, e.g. Z0=0, Zinp=0)
%       -1  - open end (Z0>>ZT, e.g. ZT=0)

% Blood properties
mu = 0.0035;    % Pa*s, dynamic viscosity
rho = 1060;     % kg/m^3, blood density
nu = mu/rho;    % m^2/s, dynamic viscosity of blood

% Terminal mean pressure
% Quereshi 2014: mean pressure = 10mmHg at artery-vein boundary set at
%   Rmin=50micron
p_terminal = 10*133.32; % Pa, or 10 mmHg

% Terminal reflection coefficient
% Hollander-2001: pulmonary circulation is of open-end type reflector
% although in McDonald book (6th ed): pulsatililty remais at capillary 
% level, 5-10% of that in arteries
GammaVal = -1;

%% Flow and Pressure Parameters
% waveform data from [Zambrano 2018]
load FlowBZ.dat;
time = FlowBZ(:,1);
flow = FlowBZ(:,2); 

% start from the 4th generation (after RIA in Olufsen 2012) of symm. tree
flow = flow./2^3; %4th gen from MPA

% define steady input flow
q_parent =  mean(flow);

% time steps
Nf = size(time,1);
T = time(Nf) - time(1);
Nt = 125;

% frequency modes, filter higher modes
NumModes = 10; 

%% =================== RUN HOMEOSTATIC OPTIMIZATION =======================
% - Optimization of Mass, Geometry and Material parameters
% - Steady solution for the symmetric tree
%   - Steady flow and mid-pressure for Mass Optimization
%   - R>Rmin, pTermSteady(Newgen) is given
%   - Mt is the mass of Collagen and SMC, Me is the mass of elastin 
%   - SMCtoCOL is ratio of content of smc to collagen
%   - Mtotal = Mc + Ms + Me; Mtotal = Mt*(5/4+SMCtoCOL)/(SMCtoCOL+1) + Me (VF?)

%-------Step 0: Optimization for initial number of generations-------------
disp('Start Homeostatic Optimization');
[Mt,Me,SMCtoCOL,Radius,Length,Table,p_mid]=...
                            TreeOptimization2(q_parent,p_terminal,N_gen); 

% output                        
disp(['    N_gen=',num2str(N_gen),': Radius=',num2str(Radius(N_gen)),...
    ', pTermSteady=',num2str(Table(N_gen,3))]);
% quick verification test
if Table(N_gen,3)~=p_terminal
    disp('ERROR: Step 0 - terminal pressure is not assigned correctly')
end

% find the number of generations related to Rmin
kk=1;
while (kk<=N_gen)&&(Radius(kk)>Rmin)
    kk=kk+1;
end
if (kk==N_gen+1) % initial number of gen. gets R>Rmin for all gen.
    Newgen=N_gen;
else
    Newgen=kk-1;
end

%--------Step 1: for new number of generations rerun optimization----------
% VF warning: may need iterative process to find number of generations
if N_gen~=Newgen
    [Mt,Me,SMCtoCOL,Radius,Length,Table,p_mid]=...
                            TreeOptimization2(q_parent,p_terminal,Newgen); 
    % output 
    disp(['Newgen=',num2str(Newgen),': Radius=',num2str(Radius(Newgen)),...
        ', pTermSteady=',num2str(Table(Newgen,3))]);
    % quick verification test
    if Table(Newgen,3)~=p_terminal
        disp('ERROR: Step 1 - terminal pressure is not assigned correctly')
    end
end
disp('End Homeostatic Optimization');

% preallocate material and geometry parameters
Thickness = zeros(1,Newgen);
YoungMod_tt = zeros(1,Newgen);
YoungMod_zz = zeros(1,Newgen);
nu_tz = zeros(1,Newgen);
nu_zt = zeros(1,Newgen);
StiffMatrix = zeros(2,2,Newgen);
ksi = zeros(1,Newgen-1);

% preallocate homeostatic values
ratio = zeros(1,Newgen);
sigma_h = zeros(1,Newgen);
shear = zeros(1,Newgen);
Pmid = zeros(1,Newgen);

% preallocate steady state solution
qSteady = zeros(1,Newgen);
pInpSteady = zeros(1,Newgen);
pTermSteady =zeros(1,Newgen);
HydRes = zeros(1,Newgen);

% assign given input mean-flow and terminal mean-pressure
qSteady(1) = q_parent; 
pTermSteady(Newgen) = p_terminal; 

%% Find unstressed radius
RzeroP = zeros(1,Newgen);
for kk=1:Newgen
    RzeroP(kk)=ZeroP(Me(kk),Mt(kk),Radius(kk),p_mid(kk),0);
end

%% update parameters and get steady solution
for k=1:Newgen
    % steady solution
    % pInpSteady(k) = pTermSteady(k)+ HydRes(k)*qSteady(k);
    qSteady(k) = Table(k,2);
    pTermSteady(k) = Table(k,3);
    HydRes(k) = Table(k,4);

    % geometry and material parameters
    % use mid-pressure
    [YoungMod_tt(k),YoungMod_zz(k),nu_tz(k),nu_zt(k),StiffMatrix(:,:,k)] =...
           YoungMod_2(Radius(k), Me(k),Mt(k),1,1,p_mid(k),SMCtoCOL(k));
       
    Thickness(k) =(1/(SMCtoCOL(k)+1)*(1.25)*Mt(k) + ...
         SMCtoCOL(k)/(SMCtoCOL(k)+1)*Mt(k) + Me(k))/(0.3*rho_w);

    % for output
    Pmid(k) = p_mid(k);
    ratio(k) = Thickness(k)/(2*Radius(k));
    sigma_h(k) = Pmid(k)/(2*ratio(k));       % homeostatic stress
    shear(k) = 4*mu*qSteady(k)/(pi*Radius(k)^3);% homeostatic shear stress
end  

% get input mean pressure
pInpSteady(1) = pTermSteady(1)+ HydRes(1)*qSteady(1);
for k=2:Newgen
    pInpSteady(k) = pTermSteady(k-1);
end

% get Murray's law exponent for output
for k=1:Newgen-1
    ksi(k) = 1/(log2(Radius(k))-log2(Radius(k+1))); 
end

%% ======================= RUN PULSATILE FLOW =============================
%   - deformable wall Womersley solution

disp('Start Pulsatile Solutions');

% get flow at root vessel in frequency domain
% QnInp(1) - steady flow, QnInp(k>1) oscillatory flow in frequency domain
% Nf-1: don't consider the last point, it's identical to the first point
[QnInp] = FlowRateFrequencyDomain(flow,NumModes,Nf-1);

% for longitudinally constrained vessel (tethered)
Ctt(1:Newgen)=StiffMatrix(1,1,:);

[pInpTime,pTermTime,qInpTime,qTermTime,cn,Alpha,zinp,zterm,zchar,InpImpedance] ...
    = WomersleySolutionSymmetricTree(Newgen,NumModes,Nt,T,rho,...
                             Radius,Length,Thickness,Ctt,GammaVal,QnInp,...
                             qSteady,pInpSteady,pTermSteady);

% real wave speed for 1st frequency (n=2)
c_R2 = zeros(1,Newgen);c_R10 = zeros(1,Newgen);
for k=1:Newgen
    c_R2(k) = 1.0/real(cn(2,k)^(-1));
    c_R10(k) = 1.0/real(cn(10,k)^(-1));
end

% Moens-Korteweg pulse wave velocity
c0_MK = zeros(1,Newgen); 
for k=1:Newgen
    c0_MK(k) = sqrt(Thickness(k)*YoungMod_tt(k)/(2*rho_w*Radius(k)));
end

% Ett pulse wave velocity
c0_tt = zeros(1,Newgen); 
for k=1:Newgen
    c0_tt(k) = sqrt(Thickness(k)*Ctt(k)/(2*rho_w*Radius(k)));
end
%% Estimate delta parameter (Vasilina's draft of CMM-Womersley paper)
% it checks if the longitudinal velocity is much smaller than the pulse
% wave velocity, delta<<1, needed to satisfy long wave approximation
% I estimated based on the maximum longitudinal velocity, computed from
% the input flow
delta1 = zeros(1,Newgen);delta9 = zeros(1,Newgen); deltaMK = zeros(1,Newgen);
for k=1:Newgen
    maxqpulse = max(qInpTime(:,k)) - qSteady(k);
    delta1(k) = maxqpulse/(pi*Radius(k)^2*c_R2(k));
    delta9(k) = maxqpulse/(pi*Radius(k)^2*c_R10(k));
    deltaMK(k) = maxqpulse/(pi*Radius(k)^2*c0_MK(k));

%         if (delta(k) >= 0.04)
%     if (delta(k) >= 0.1)
%         disp(['WAWRNING: not small delta - Womersley solution is not valid in ',num2str(k),' generation!']);
%     end
end

%% Get Root Impedance
[Zn,z] = GetImpedance (qInpTime(:,1),pInpTime(:,1),T,Nf,Nt,NumModes);


disp('End Pulsatile Solutions');


% total for the root vessel zinp,zterm,zchar,InpImpedance
zinp_total=pInpSteady(1)./qSteady(1) + zinp(:,1);
zterm_total=pTermSteady(1)./qSteady(1) + zterm(:,1);

%% ========================== FIGURES =====================================
% Output for the presentatoins and paper

%% ************************** homeostatic optimization results ************

% mass fraction
figure(1);
hold on;
for i=1:Newgen
    [e, c, m] = mass_fracs(2*Radius(i));
    semilogx(2*Radius(i),e,'ob');
    semilogx(2*Radius(i),m,'or');
    semilogx(2*Radius(i),1 - e - m,'og');
    grid on;
end
xlabel('log(D[\mum])'); ylabel('Mass fraction');
  legend('elastin','SMC','collagen');
  hold off;


% varable for generation number 
k = 1:1:Newgen;  

h1 = figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
%     set(gca,'xscale');
set(gca,'XDir','reverse');
xlabel('D (cm)');ylabel('E_{\theta\theta} (kPa)');hold on;
ylim([30 60]);
plot(2*Radius(k)*100,YoungMod_tt(k)/1000,'o-','LineWidth',2); %,'LineWidth',1.5);
FigModify(h1,'E')

h2 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'XDir','reverse');xlabel('D (cm)');ylabel('\sigma_h (kPa)');hold on;
ylim([5 15]);
plot(2*Radius(k)*100,sigma_h(k)/1000,'o-','LineWidth',2); %,'LineWidth',1.5);
box on;
print(h2, '-dpdf', ['sigma_h_var','.pdf']);
% FigModify(h2,'sigma_h')

h3 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'XDir','reverse');xlabel('D (cm)');ylabel('\tau (Pa)');hold on;
ylim([1.00 1.30]);
plot(2*Const.Radius(k)*100,Const.shear(k),'o-','LineWidth',2); %,'LineWidth',1.5);
plot(2*Radius(k)*100,shear(k),'o-','LineWidth',2); %,'LineWidth',1.5);
box on;
legend('Case 1','Case2');
print(h3, '-dpdf', ['shear_var','.pdf']);
% FigModify(h3,'shear')


h4 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'XDir','reverse');xlabel('D (cm)');ylabel('p_{mid} (mmHg)'); hold on;
ylim([8 13]);
plot(2*Radius(k)*100,Pmid(k)/133.32,'o-','LineWidth',2); %,'LineWidth',1.5);
% FigModify(h4,'Pmid')
box on;
print(h4, '-dpdf', ['Pmid_var','.pdf']);

h5 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');grid on;
plot(2*Radius(k)*100,ratio(k),'o-','LineWidth',2); hold on;
errorbar(0.002136*100,0.162/2.136,-0.027/2.136,0.027/2.136,'d','LineWidth',2);
errorbar(0.000250*100,0.0815,-0.0035,+0.0035,'dk','LineWidth',2);
set(gca,'XDir','reverse');xlabel({'D (cm)'},'FontSize',12);ylabel({'h/D'},'FontSize',12);hold on;
legend({'Thickness-to-diameter ratio','[Li-2012]','[Rol-2017]'},'FontSize',12,'Location','SouthWest');
ylim([0.00 0.12]);
box on;
print(h5, '-dpdf', ['ratio_var','.pdf']);
% FigModify(h5,'ratio')

h6 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'xscale');xlabel('gen.');ylabel(' Radius exponent \xi');hold on;
ylim([2.5 3.1]);
plot(1:Newgen-1,Const.ksi(1:Newgen-1),'o-','LineWidth',2); %,'LineWidth',1.5);
plot(1:Newgen-1,ksi(1:Newgen-1),'o-','LineWidth',2); %,'LineWidth',1.5);
xlim([0 Newgen]);
legend('Case 1','Case 2');
% FigModify(h6,'xi')
box on;
print(h6, '-dpdf', ['xi_var','.pdf']);
% suptitle(['No. generations = ',num2str(Newgen),' for D_{min} = ',num2str(2*Rmin*100),'cm']);

% hold off;

%% Postprocessing geometrical values vs generations

% figure(11);

% plot tree parameters R,L,H,E versus generations
k = 1:1:Newgen; 

h13 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');hold on;grid on;
plot(k, 2*Radius(k).*100,'o-','LineWidth',2); 
plot(1 ,2*R0*100,'*','LineWidth',2);
plot(1 ,0.271,'*','LineWidth',2);
plot(1 ,0.416,'*','LineWidth',2);
plot(1 ,0.734,'*','LineWidth',2);
legend('Optimization results','3rd gen. [Olufsen-2012]','12th order [Huang-1996]','13th order [Huang-1996]','14th order [Huang-1996]');
xlabel('gen.');ylabel('D (cm)');
% FigModify(h13,'D_to_gen');
box on;
print(h13, '-dpdf', ['D_to_gen_var','.pdf']);


hLength = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');hold on;grid on;
plot(k, Length(k).*100,'o-','LineWidth',2); 
plot(1 ,2*R0*100,'*','LineWidth',2);
plot(1 ,0.271,'*','LineWidth',2);
plot(1 ,0.416,'*','LineWidth',2);
plot(1 ,0.734,'*','LineWidth',2);
legend('Optimization results','3rd gen. [Olufsen-2012]','12th order [Huang-1996]','13th order [Huang-1996]','14th order [Huang-1996]');
xlabel('gen.');ylabel('D (cm)');
% FigModify(h13,'D_to_gen');
box on;
print(hLength, '-dpdf', ['Length_var','.pdf']);

% subplot(2,2,2); set(gca,'xscale');  xlim([0 Newgen]); hold on;
% plot(k, Thickness(k).*100,'o-','LineWidth',2); 
% plot(1 ,H0*100,'*','LineWidth',2);
% legend('homeostatic values','H_0=0.1*R_0 (Banks)');
% xlabel('gen.');ylabel('h (cm)');
% 
% subplot(2,2,3); set(gca,'xscale'); xlim([0 Newgen]); hold on;
% plot(k, Length(k).*100,'o-','LineWidth',2);   % factor 5 from Olufsen
% plot(k, Length(k).*2*100,'--','LineWidth',2); %back to Olufsen curve
% plot(1, L0*100,'*','LineWidth',2); 
% plot(1, 1.807,'*','LineWidth',2); 
% % legend('homeostatic values','L_0 (RIA, Olufsen)','L (12 order, Huang)');
% legend('homeostatic values=L_{Oluf}/2','L_{Oluf}(r)-Olufsen fitting','L_0 (RIA, Olufsen)','L (12 order, Huang)');
% xlabel('gen.');ylabel('l (cm)');
% 
h14 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  6 5],'PaperPositionMode','auto');hold on;grid on;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
ylim([1.1 1.4]);
eta01=1.2; eta02=1.3;
eta=ones(1,Newgen); %(A1+A2)/A0=2A1/A0=2R1^2/R0^2
for m=2:Newgen
    eta(m)=2*Radius(m)^2/Radius(m-1)^2;
end
etaR1(1:Newgen)=eta01;etaR2(1:Newgen)=eta02;
plot(1:Newgen-1, eta(2:Newgen),'o--','LineWidth',2); 
plot(1:Newgen-1,etaR1(2:Newgen),'r--','LineWidth',2);
plot(1:Newgen-1,etaR2(2:Newgen),'r--','LineWidth',2);
legend('homeostatic values','open-end reflection bounds');
xlabel('gen.');ylabel('Daughter-to-parent area ratio, \eta');
% FigModify(h14,'eta');
box on;
print(h14, '-dpdf', ['eta_var','.pdf']);
% 
% suptitle(['No. generations = ',num2str(Newgen),' for D_{min} = ',num2str(2*Rmin*100),'cm']);
% 
% hold off;

%% Postprocessing properties vs generations

% plot tree parameters R,L,H,E versus generations
PWV1(1:Newgen)=PwvB; PWV2(1:Newgen)=PwvM;

h7 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');grid on;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k, Const.YoungMod_tt(k)./1000,'o-','LineWidth',2,'Color', [0 114/255 189/255]); 
plot(k, YoungMod_tt(k)./1000,'o-','LineWidth',2,'Color', [217/255, 83/255, 25/255]); 
plot(k, Const.YoungMod_zz(k)./1000,'+-','LineWidth',2,'Color', [0 114/255 189/255]); 
plot(k, YoungMod_zz(k)./1000,'+-','LineWidth',2,'Color', [217/255, 83/255, 25/255]); 
legend('Case 1: E_{\theta\theta}',' Case 2: E_{\theta\theta}','Case 1: E_{zz}','Case 2: E_{zz}','Location','SouthWest');
xlabel('gen.');ylabel('Young modulus (kPa)');
% FigModify(h7,'YoungMod')
box on;
print(h7, '-dpdf', ['YoungMod_var','.pdf']);


h8 = figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
nu0(1:Newgen)=0.5;
ylim([0.1 0.7]);
plot(k, nu_tz(k),'o-','LineWidth',2); 
plot(k, nu_zt(k),'o-','LineWidth',2); 
plot(k, nu0,'--','LineWidth',2); 
legend('\nu_{{\theta}z}','\nu_{z\theta}','\nu=0.5','Location','NorthEast');
xlabel('gen.');ylabel('Poisson ratio');
FigModify(h8,'PoissonRatio')


% h9 = figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');
% set(gca,'xscale'); xlim([0 Newgen]); hold on;
% plot(k, Thickness(k).*YoungMod_tt(k)./1000,'o--','LineWidth',2); 
% plot(1, H0*E0/1000,'*','LineWidth',2); 
% plot(1, EhZ/1000,'*','LineWidth',2); 
% legend('homeostatic values','? E_{\theta\theta}(est.)*H0','? Eh (Zambrano)','Location','NorthEast');
% xlabel('gen.');ylabel('E_{\theta\theta}*h (kPa*m)');
% FigModify(h9,'Eh')
%%%%%%%%%% Yen et al's data
% EhR_Yen;
%%%%%%%%%%

h10 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');grid on;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
% ylim([4 12]);
EhrKall(1:Newgen)=EhrK; EhrYall(1:Newgen)=EhrY; EhrQall(1:Newgen)=EhrQ;
plot(k, Const.Thickness(k).*Const.YoungMod_tt(k)./Const.RzeroP(k)./1000,'o--','LineWidth',2); 
plot(k, Thickness(k).*YoungMod_tt(k)./RzeroP(k)./1000,'o--','LineWidth',2); 
plot(k, EhrKall./1000,'--','LineWidth',2); 
% plot(k, EhrQall./1000,'--');
plot(k, EhrYall./1000,'-.','LineWidth',2); 
errorbar(1,EhrEx_mean/1000,abs(min(EhrEx_err/1000)),abs(max(EhrEx_err/1000)),'o','LineWidth',2); 
legend({'Case 1','Case 2','Eh/R0 [Krenz-2003]','Eh/R0 [Yen-1990]','E_{\theta\theta}h/R0 (MSU Experiments)'},'Location','SouthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('E_{\theta\theta}h/R0 (kPa)','FontSize',12);box on;
axis([0 19 0 12])
box on;
print(h10, '-dpdf', ['Ehr_var','.pdf']);

k = 1:1:Newgen; 
h11 = figure;hold on; grid on;%figure('Units','inches','Position',[4 4  7 6],'PaperPositionMode','auto');grid on;
set(gca,'xscale');  xlim([0 Newgen]); hold on;
plot(k, c_R2,'o-','LineWidth',2);
plot(k, c_R10,'o-','LineWidth',2); 
% plot(k ,c0_MK,'.-r','LineWidth',2); 
plot(k ,PWV1,'--','LineWidth',2);
plot(k ,PWV2,'-.','LineWidth',2);
% plot(1 ,PwvB,'*','LineWidth',2); 
% plot(1 ,PwvM,'*','LineWidth',2); 
legend({'at 1st mode \omega=2\pi/T','at 9th mode \omega=2\pi9/T',...
 'Wave velocity [Banks-1978]',...
 'Wave velocity [Milnor-1969]'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('Pulse wave velocity (m/s)','FontSize',12);
FigModify(h11,'PWV_var')

h12 = figure('Units','inches','Position',[4 4  4.6 3.8],'PaperPositionMode','auto');grid on;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
ub(1:Newgen)=1;tb(1:Newgen)=0.1;
plot(k,delta1(k),'o-','LineWidth',2); 
plot(k,deltaMK(k),'bs-','LineWidth',2); 
% plot(k,ub(k),'--');
% plot(k,tb(k),'-.');
legend('delta=max(longVel)/c_{r1}','deltaMK=max(longVel)/c_{MK}','upper bound','treshold: delta<<1','Location','NorthEast');
legend('delta=max(longVel)/c_{r1}','deltaMK=max(longVel)/c_{MK}','Location','NorthEast');
xlabel('gen.');ylabel('delta');
FigModify(h12,'delta')
hold off;
% 

%% mass fraction
h701 = figure('Units','inches','Position',[4 4  6 5],'PaperPositionMode','auto');
DD = logspace(-2.1,4,100);
[e, c, m] = mass_fracs(DD/100);
semilogx(DD,e,'b','linewidth',2);hold on;
semilogx(DD,m,'r','linewidth',2);
semilogx(DD,1 - e - m,'Color',[.61 .11 .9],'linewidth',2);
semilogx([3200/10^4 3200/10^4],[0 1],'k','linewidth',2);
semilogx([2000/10^4 2000/10^4],[0 1],'k','linewidth',2);
% set(h701,'linewidth',2);
grid on;
set(gca,'XDir','reverse');
ylim([0 1.0]);
xlim([0.008 1.0]);
xlabel('D (cm)'); ylabel('Mass fraction');
legend('elastin','SMC','collagen','Location','NorthEast');
hold off;
FigModify(h701,'mass_fractions')


%% ************************** pulsatile hemodynamics results ************
%------------------------------PWV-----------------------------------------
hp1=figure;
set(gca,'xscale');  xlim([0 Newgen]); hold on;
plot(k, c_R2,'o-','LineWidth',2); 
plot(k, c_R10,'o-','LineWidth',2); 
plot(k ,c0_MK,'.-k','LineWidth',2);
% plot(k ,c0_tt,'--k','LineWidth',2);
plot(1 ,PwvB,'*','LineWidth',2);
plot(1 ,PwvM,'*','LineWidth',2);
legend({'at 1st mode \omega=2\pi/T','at 9th mode \omega=9\times2\pi/T',...
 'c^{MK} with E_{\theta\theta}','[Banks-1978]',...
 '[Milnor-1969]'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('Pulse wave velocity (m/s)','FontSize',12);
grid on;
axis([0 Newgen 0 3]);
box on;
print(hp1, '-dpdf', ['PWV_var','.pdf']);
% FigModify(hp1,'PWV')

%------------------------------Delta---------------------------------------
hp2=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k,delta1(k),'o-','LineWidth',2);
plot(k,delta9(k),'o-','LineWidth',2);
% legend({'delta=max(longVel)/c_{r1}','deltaMK=max(longVel)/c_{MK}'},'Location','NorthEast','FontSize',12);
legend({'{\delta} at 1st mode','{\delta} at 9th mode'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('{\delta}','FontSize',12);
grid on;
% FigModify(hp2,'delta')
box on;
print(hp2, '-dpdf', ['delta_var','.pdf']);

%--------------------------Womersley number--------------------------------
hp3=figure;
set(gca,'xscale'); xlim([0 Newgen]); hold on;
plot(k,Alpha(2,k),'o-','LineWidth',2);
plot(k,Alpha(10,k),'o-','LineWidth',2);
legend({'{\alpha} at 1st mode  \omega=2\pi/T','{\alpha} at 9th mode \omega=9\times2\pi/T'},'Location','NorthEast','FontSize',12);
xlabel('gen.','FontSize',12);ylabel('Womersley number','FontSize',12);
grid on;
box on;
print(hp3, '-dpdf', ['WomresleyNo_var','.pdf']);
% FigModify(hp3,'WomresleyNo')

% time points
t = linspace(0,T-T/Nt,Nt);
%--------------------------Flow and Pressure-------------------------------
hp4=figure;
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
xlabel('Time (s)'); ylabel('Flow (cm^3/s)'); grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T -5 50]);
box on
% set(gca,'FontSize',16)
% set(gca,...
% 'Units','normalized',...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',16,...
% 'FontName','Calibri')
% print(hp7, '-dpdf', [name,'.pdf']);
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
xlabel('Time (s)');  ylabel('Presure (mmHg)');   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T 5 25]);
% box on
% set(gca,'FontSize',16)
% set(gca,...
% 'Units','normalized',...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',16,...
% 'FontName','Calibri')

% print(hp4, '-dpdf', ['PressureFlow_var','.pdf']);
box on;
print(hp4, '-dpdf', ['PressureFlow_var','.pdf']);

%------------------------ Root pressure------------------------------------
hp5=figure;
plot(t,pInpTime(:,1)/133.32,'-','LineWidth',2);
hold on
plot(t,pTermTime(:,1)/133.32,'-','LineWidth',2);
legend({'Input p at 1st gen.','Terminal p at 1st gen.'},'FontSize',12);
xlabel('Time (s)','FontSize',12);  ylabel('Presure (mmHg)','FontSize',12);   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T 5 25]);
box on;
print(hp5, '-dpdf', ['pTerm1_var','.pdf']);
%------------------------ Root flow----------------------------------------
hp6=figure;
plot(t,qInpTime(:,1)*10^6,'-','LineWidth',2);
hold on
plot(t,qTermTime(:,1)*10^6,'-','LineWidth',2);
legend({'Input q at 1st gen.','Terminal q at 1st gen.'},'FontSize',12);
xlabel('Time (s)','FontSize',12);  ylabel('Flow (cm^3/s)','FontSize',12);   grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
axis([0 T -5 50]);
box on;
print(hp6, '-dpdf', ['qInp1_var','.pdf']);

%---------------------- Impedance in time domain---------------------------
hp8=figure;
% plot(t,z./(133.32*10^6),'--','LineWidth',2);  hold on
plot(t,zinp_total./(133.32*10^6),'-','LineWidth',2);  hold on
plot(t,zterm_total./(133.32*10^6),'-.','LineWidth',1);  hold on
plot(t,zchar(:,1)./(133.32*10^6),'--','LineWidth',1);  hold off
grid on;
set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
legend({'Z_{inp}','Z_{term}','Z_c-(pulsatile)'},'FontSize',12);
xlabel('Time (s)','FontSize',12); ylabel('Impedance(mmHg*s/cm^3)','FontSize',12);
axis([0 T -4 10]);
box on;
print(hp8, '-dpdf', ['ImpTime_var','.pdf']);
%---------------------- Impedance in frequency domain----------------------
hp7=figure;
kk = 0:1:NumModes-1;
% suptitle('Input impedance in frequency domain');
subplot(1,2,2); plot(kk,angle(Zn),'o','LineWidth',2); grid on;
box on
% set(gca,'FontSize',16)
% set(gca,...
% 'Units','normalized',...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',16,...
% 'FontName','Calibri')
% set(gca,'ytick',-pi:pi/4:pi); 
% set(gca,'yticklabel',{'-\pi','-\pi/2','-\pi/4','0','\pi/4','\pi/2','\pi'}); 
xlabel('Frequency modes'); ylabel('Phase Angle');
subplot(1,2,1); plot(kk,abs(Zn)./(133.32*10^6),'o','LineWidth',2); grid on;
set(gca,'xtick'); xlabel('Frequency modes'); ylabel('Modulus (mmHg*s/cm^3)');
% FigModify(hp7,'ImpFreq')
box on;
print(hp7, '-dpdf', ['ImpFreq_var','.pdf']);
