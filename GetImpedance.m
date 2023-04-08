%--------------------------------------------------------------------------
% Given flow and pressure waveforms, compute the impedance function.
% - Input arrays q(t) and p(t) should be preprocessed:
%   - have the same size, Nf,(interpolate if needed);
%   - adjust data to the same time period, T;
%   - synchronize in time (from ECG data) to have physiological phase-lag.
%   - have the same system of units:
%       - in centimetre–gram–second, q[cm^3/s] and p[dyn/cm^2];
%       - in IS, q[m^3/s] and p[Pa], note 1Pa=10dyn/cm^2;
% - Here p(t),q(t),z(t) in time-domain; Pn(k),Qn(k),Zn(k) in freq-domain
%   - use FFT and IFFT.
% - The script filters q and p to first several modes (about 10)
%   - make sure that q has no near-zero k-th mode values, Qn(k); if so,
%   filter them out (needed to ensure no-division-by-zero in Zn(k)).
% - The script tests the resulted impedance:
%   - convolution integral of q(t)z(t-t1) should give p(t);
% - The number of  modes, NumModes, controls the accuracy of filtering
%   - The higher NumModes the higher frequencies are in z(t), which may 
%   bring concerns to flowsolver convergency and stability issues. Thus  
%   NumModes should be optimal: low enough to have smoother z(t) but not 
%   to loose p,q accuracy;
%   - The script has reconstructed p&q to reference the filtering accuracy.
%
% (by Vasilina, 07-13-2018)
%--------------------------------------------------------------------------
% %% Prepare workspace
% clc; clear all; close all;
% format long; workspace
% 
% %% Inputs and Parameters
% 
% % Read input files
% % q[cm^3/s] and p[dyn/cm^2]
% load flow.inp;
% load pressure.inp;
% q = flow(:,2);
% p = pressure(:,2);
% t = flow(:,1); % = pressure(:,1)
% Nf = size(t,1); 
% T = t(Nf)-t(1);
% 
% % Parameters
% NumModes = 13; %number of modes in filtering
% Nt = Nf-1; % time increments for z(t) output, should have Nt > 2*NumModes-2
% epsilon = 1.e-4; %tollerance for Qn mode trancation

function [Zn,z] = GetImpedance (q,p,T,Nf,Nt,NumModes)
epsilon = 1.e-8; %tollerance for Qn mode trancation

%% Computations

% FFT: from time to frequency domain  
% Qn(1),Pn(1)-steady, Qn(k>1),P(k>1)-oscillatory
% Nf-1: don't consider the last point, it's identical to the first point
[Qn] = FilterAndFrequencyDomain(q,NumModes,Nf-1);
[Pn] = FilterAndFrequencyDomain(p,NumModes,Nf-1);

% Check mode truncation
for k=1:NumModes
    if (abs(Qn(k))<epsilon)
        disp(['WARNING: too small Qn at k=',num2str(k),'; make NumModes<k.']);
        break;
    end
end

% Get impedance Zn(1:NumModes) in frequency domain
% already truncated (as Q and P are filtered)
Zn = Pn./Qn;

% Get impedance z(1:Nt) in time domain
% Note that output is real(z)
z = FilterAndTimeDomain(NumModes,Nt,Zn);

%Test with convolution integral
pconv(1:Nt) = cconv(q(1:Nt),z,Nt)./Nt;

%Test p,q reconstraction after filtering
prec = FilterAndTimeDomain(NumModes,Nt,Pn);
qrec = FilterAndTimeDomain(NumModes,Nt,Qn);


% %% Output and Figures
% tt = linspace(0,T-T/Nt,Nt);
% kk = 1:1:NumModes;
% 
% % Output z(t)
% % Note: for OldFlowsolevr z should be scaled by factor Nt: z=z*Nt
% fid = fopen('Impedance.dat','w');
% y = zeros(2,Nt+1);
% y(1,1:Nt) = tt(1,:); y(1,Nt+1) = T;
% y(2,1:Nt) = z;      y(2,Nt+1) = z(1);
% fprintf(fid,'%10.6f %12.8f\n',y);
% fclose(fid);
% 
% % Plot data and results
% 
% figure(600)
% plot(tt,q(1:Nt)*10^6,'r-'); hold on;
% plot(tt,qrec*10^6,'k--');  hold on;
% hold off;
% grid on;
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% title('Input waveforms: Flow(cm^3/s)');
% xlabel('Time (s)'); ylabel('Waveforms');
% legend('Original Flow','Flow-filtered');
% 
% figure(601)
% plot(tt,p(1:Nt)./133.32,'b-'); hold on;
% plot(tt,prec./133.32,'g--');
% hold off;
% grid on;
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% title('Input waveforms: Pressure(mmHg)');
% xlabel('Time (s)'); ylabel('Waveforms');
% legend('Original Pressure','Pressure-filtered');
% 
% figure(602);
% suptitle('Flow rate in frequency domain');
% subplot(1,2,1); plot(kk,angle(Qn),'o'); grid on;
% set(gca,'ytick',-pi:pi/2:pi); 
% set(gca,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}); 
% xlabel('Frequency modes'); ylabel('Phase Angle');
% subplot(1,2,2); plot(kk,abs(Qn),'o'); grid on;
% set(gca,'xtick'); xlabel('Frequency modes'); ylabel('Modulus');
% 
% figure(603);
% plot(tt,z./(133.32*10^6),'-','LineWidth',2);  
% grid on;
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% legend({'Impedance in time domain'},'FontSize',12);
% xlabel('Time (s)','FontSize',12); ylabel('Impedance(mmHg*s/cm^3)','FontSize',12);
% axis([0 T -4 10]);
% 
% figure(604);
% suptitle('Impedance in frequency domain');
% subplot(1,2,1); plot(kk,angle(Zn),'o','LineWidth',2); grid on;
% set(gca,'ytick',-pi:pi/2:pi); 
% set(gca,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}); 
% xlabel('Frequency modes'); ylabel('Phase Angle');
% subplot(1,2,2); plot(kk,abs(Zn),'o','LineWidth',2); grid on;
% set(gca,'xtick'); xlabel('Frequency modes','FontSize',12); ylabel('Modulus','FontSize',12);
% 
% figure(605);
% plot(tt,p(1:Nt),'b-');
% hold on;
% plot(tt,pconv,'m--');
% hold off;
% grid on;
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% title('Test pressure');
% xlabel('Time (s)'); ylabel('Pressure (dyn/cm^2)');
% legend('Inpute pressure','Convoluted pressure');
