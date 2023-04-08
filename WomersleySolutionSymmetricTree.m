% updated 07-08-18
%   - added tethered solution
%   - corrected for Ctttt instead of E,sigma
%       C = [Ctttt,Cttzz;Cttzz,Czzzz]./((Mt+Me)/(0.3*rho_w))

function [pInpTime,pTermTime,qInpTime,qTermTime,WaveSpeed,Alpha,zinp,zterm,zchar,InpImpedance] ...
    = WomersleySolutionSymmetricTree(Ngen,NumModes,Nt,T,rho,...
                            Radius,Length,Thickness,Stiff,GammaVal,...
                            QnInp,qSteady,pInpSteady,pTermSteady)

global rho_w mu
% output input impedance, pressure and flow at each generation
% symmetric tree simplification
% Radius(k) = R0*k

% global parameters:
% mu - blood kinematic viscosity;


%--------- Preallocation and initialization-------------
% Impedance,WaveSpeed are zero at first mode

% frequency domain
GammaList = zeros(NumModes,Ngen);
CharactImpedance= zeros(NumModes,Ngen);
InpImpedance = zeros(NumModes,Ngen);
TermImpedance = zeros(NumModes,Ngen);
WaveSpeed = zeros(NumModes,Ngen);
omegan = zeros(NumModes,1);
GammaT = zeros(1,NumModes);
ZTerm = zeros(1,NumModes);
Alpha = zeros(NumModes,Ngen);

% time domain
zinp = zeros(Nt,Ngen);
zterm = zeros(Nt,Ngen);
zchar = zeros(Nt,Ngen);

%==========================================================================
%---------------- Backward Computation of Impedance -------------
%==========================================================================
%--------------- Set Frequency Modes ---------------
% omegan(1)=0,  omegan(2)=2pi/T,...omegan(n)=2pi(n-1)/T
for k=1:NumModes
	omegan(k) = 2.0*pi*(k-1)/T;
end

% terminal reflection coef for all modes (first is not involved ?)
GammaT(1:NumModes) = GammaVal;
%-------------------- initialization for last generation -----------------
% reflection coefficient
GammaList(:,Ngen) = GammaT(:);

%get Womersley coefficients
% [cn,Mn,gn,c_Rn] = WomersleySolutionCoeff(...
%                   NumModes,omegan,mu,rho,...
%                   Radius(Ngen),YoungMod(Ngen),sigma,Thickness(Ngen),rho_w);
 
%tethered wall (longitudinally constrained)      
[cn,Mn,gn,c_Rn,alphan] = WomersleySolutionCoeffTethered(...
                  NumModes,omegan,mu,rho,...
                  Radius(Ngen),Stiff(Ngen),Thickness(Ngen));
WaveSpeed(:,Ngen) = cn;
Alpha(:,Ngen) = alphan;

%get characteristic impedance
[Z0] = CharacteristicImpedanceOneVessel(NumModes,Radius(Ngen),rho,cn,Mn,gn); 
CharactImpedance(1,Ngen) = 0.0;   
CharactImpedance(2:NumModes,Ngen) = Z0(2:NumModes); 

%get input impedance
[Zinp] = InputImpedanceOneVessel(NumModes,omegan,Length(Ngen),Z0,cn,GammaT); 
InpImpedance(1,Ngen) = 0.0; 
InpImpedance(2:NumModes,Ngen) = Zinp(2:NumModes);         

%get terminal impedance
[ZT] = TerminalImpedanceVessel(NumModes,Z0,GammaT);
TermImpedance(1,Ngen) = 0.0;  
TermImpedance(2:NumModes,Ngen) = ZT(2:NumModes);         

%get time domain functions
[z0total,z0oscil] = FilterAndIFFT(NumModes,Nt,CharactImpedance(:,Ngen));
[zinpTotal,zinpOscil] = FilterAndIFFT(NumModes,Nt,InpImpedance(:,Ngen));                                   
[ztermTotal,ztermOscil] = FilterAndIFFT(NumModes,Nt,TermImpedance(:,Ngen)); 
zinp(:,Ngen) = zinpTotal(:);
zterm(:,Ngen)= ztermTotal(:);  
zchar(:,Ngen)= z0total(:);

% -------Backward loop for the rest of generations----------
% all segments are identical, solve just for one segment per generation
for k = Ngen-1:-1:1
    % compute symmetrical terminal impedance from (k+1)th generation
%     ZT = zeros(1,NumModes);
%     ZT(2:NumModes) = InpImpedance(2:NumModes,k+1)/2.0;
    [ZT] =  TerminalImpedanceJoint(NumModes,InpImpedance(:,k+1),InpImpedance(:,k+1));   
    TermImpedance(1,k) = 0.0; 
    TermImpedance(2:NumModes,k) = ZT(2:NumModes);

%     %get Womersley coefficients
%     [cn,Mn,gn,c_Rn] = WomersleySolutionCoeff(...
%         NumModes,omegan,mu,rho,...
%         Radius(k),YoungMod(k),sigma,Thickness(k),rho_w);
    
    %tethered wall (longitudinally constrained)      
    [cn,Mn,gn,c_Rn,alphan] = WomersleySolutionCoeffTethered(...
                  NumModes,omegan,mu,rho,...
        Radius(k),Stiff(k),Thickness(k));
    
    WaveSpeed(:,k)=cn;
    Alpha(:,k)=alphan;
    
    %get Z0
    [Z0] = CharacteristicImpedanceOneVessel(NumModes,Radius(k),rho,cn,Mn,gn);
    CharactImpedance(1,k) = 0.0; 
    CharactImpedance(2:NumModes,k) = Z0(2:NumModes);

    %get reflection coefficient
    [Gamma] = ReflectionCoeff(NumModes,ZT,Z0);
    GammaList(:,k) = Gamma(:);

    %get input impedance
    [Zinp] = InputImpedanceOneVessel(NumModes,omegan,Length(k),Z0,cn,Gamma); 
    InpImpedance(1,k) = 0.0;                    
    InpImpedance(2:NumModes,k) = Zinp(2:NumModes);

    %get time domain functions
    [z0total,z0oscil] = FilterAndIFFT(NumModes,Nt,CharactImpedance(:,k));
    [zinpTotal,zinpOscil] = FilterAndIFFT(NumModes,Nt,InpImpedance(:,k));                                 
    [ztermTotal,ztermOscil] = FilterAndIFFT(NumModes,Nt,TermImpedance(:,k)); 
    zinp(:,k) = zinpTotal(:);
    zterm(:,k)= ztermTotal(:);
    zchar(:,k)= z0total(:);
end

%==========================================================================
%----------- Forward Computation of Total Pressure and Flow -----------
%==========================================================================
% hemodynamics in all segments within one generation is the same for
% symmetric tree

PnInp = zeros(NumModes,Ngen);
pInpTime = zeros(Nt,Ngen);
PnTerm = zeros(NumModes,Ngen);
pTermTime = zeros(Nt,Ngen);

% qInpRoot = zeros(Nt,1);
qInpTime = zeros(Nt,Ngen);
qTermTime = zeros(Nt,Ngen);

%-----------Initialization at root --------
% get input pressure
PnInp(1,1) = pInpSteady(1);
for i=2:NumModes
    PnInp(i,1) = InpImpedance(i,1)*QnInp(i);
end
% 1) input pressure time
[RootPressInp,pt] = FilterAndIFFT(NumModes,Nt,PnInp(:,1)); 
pInpTime(:,1) = RootPressInp;

%get Hf coefficient, for oscillations only
[Hfn] = HforwardCoeff(NumModes,omegan,Length(1),PnInp(:,1),WaveSpeed(:,1),GammaList(:,1));

% % 2) input pressure time
% x=0.0;
% [Pn]= Precovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,1,1),GammaList(:,1,1));
% Pn(1) =  pInpSteady(1);    
% [RootPressInp2,pt] = FilterAndIFFT(NumModes,Nt,Pn); 
% pInpTime(:,1,1) = RootPressInp2;

% --------- flow----------
%get root input flow 
x = 0.0;
[Qn]=Qrecovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,1),...
                CharactImpedance(:,1),GammaList(:,1)); % Ok, the same as QnInp
Qn(1) = qSteady(1); 
[qtotal,qoscil] = FilterAndIFFT(NumModes,Nt,Qn);
qInpTime(:,1) = qtotal;

%get root output flow 
x = Length(1);
[Qn]=Qrecovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,1),...
                CharactImpedance(:,1),GammaList(:,1));
Qn(1) = qSteady(1); 
[qtotal2,qoscil2] = FilterAndIFFT(NumModes,Nt,Qn);
qTermTime(:,1) = qtotal2;  % ?
         
% -------- pressure -------
% %get input pressure
% x = 0.0;
% [Pn] = Precovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,1),GammaList(:,1));
% Pn(1) = pInpSteady(1); 
% PnInp(:,1) = Pn(:);  
% [RootPressInp,pt] = FilterAndIFFT(NumModes,Nt,PnInp(:,1)); 
% pInpTime(:,1) = RootPressInp; % laready obtained

%get root terminal pressure 
x = Length(1);
[PnT] = Precovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,1),GammaList(:,1));
PnT(1) = pTermSteady(1); %PnTerm(1,1) = pTermSteady(1);
PnTerm(:,1)= PnT;   
[PTtotal,pTt] = FilterAndIFFT(NumModes,Nt,PnTerm(:,1)); 
pTermTime(:,1) = PTtotal;

%---------- Get p, q for the rest of generations-------
% in frequency domain
for k = 2:Ngen
    %get input pressure as parent terminal pressure
    PnInp(:,k) = PnTerm(:,k-1); % pressure conservation
    %pInpTime(:,k) = pTermTime(:,k-1);  
    % or
    [PTtotal,pTt] = FilterAndIFFT(NumModes,Nt,PnInp(:,k)); 
    pInpTime(:,k) = PTtotal;
            
    %get Hf coefficient, for oscillations only
    [Hfn] = HforwardCoeff(NumModes,omegan,Length(k),PnInp(:,k),WaveSpeed(:,k),GammaList(:,k));              

    %get terminal pressure
    x = Length(k);
    [PnT]=Precovery(NumModes,omegan,Length(k),x,Hfn,WaveSpeed(:,k),GammaList(:,k));
    PnT(1) = pTermSteady(k);  %PnTerm(1,k) = pTermSteady(k); 
    PnTerm(:,k) = PnT; %oscillatory
    [PTtotal,pTt] = FilterAndIFFT(NumModes,Nt,PnTerm(:,k)); 
    pTermTime(:,k) = PTtotal;
    
    %get input flow
    x = 0.0;
    [Qn]=Qrecovery(NumModes,omegan,Length(1),x,Hfn,WaveSpeed(:,k),...
                CharactImpedance(:,k),GammaList(:,k)); 
    Qn(1) = qSteady(k); 
    [qtotal,qoscil] = FilterAndIFFT(NumModes,Nt,Qn);
    qInpTime(:,k) = qtotal;
        
    % get terminal flow
    x = Length(k);
    [Qn]=Qrecovery(NumModes,omegan,Length(k),x,Hfn,WaveSpeed(:,k),...
                CharactImpedance(:,k),GammaList(:,k));
    Qn(1) = qSteady(k); 
    [qtotal2,qoscil2] = FilterAndIFFT(NumModes,Nt,Qn);
    qTermTime(:,k) = qtotal2; 
end
