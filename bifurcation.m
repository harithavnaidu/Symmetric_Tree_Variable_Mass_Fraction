function [Mt1, Mt2, lhat1, lhat2,Q_parent,P1,P2] = bifurcation(Q_parent,p_terminal,Me)

global R0 mu
%% Bifurcation model!
for i = 1:50
    Q_parent(i) = (i+20)*10^(-6)/60;
    
    err = 10;
    lhat1(i) = 1.0;
    lhat2(i) = 1.0;
    
    length_ratio = 42;
    L1 = length_ratio*R0*lhat1(i);
    L2 = length_ratio*R0*lhat2(i);
    q1 = Q_parent(i);
    q2 = Q_parent(i)/2;
    
    p_mid2 = p_terminal + 8*mu*(L2/2)*q2/(pi*(lhat2(i)*R0)^4);
    p_bifur = p_terminal + 8*mu*(L2)*q2/(pi*(lhat2(i)*R0)^4);
    p_mid1 = p_bifur + 8*mu*(L1/2)*q1/(pi*(lhat1(i)*R0)^4);
    
    Y0=[p_mid1,p_bifur,p_mid2];
    c=0;
    
    while err > 1e-6
        
        
        [Mt1(i), lhat1(i)] = mass_optimiz(p_mid1,q1,Me(1));
        
        [Mt2(i), lhat2(i)] = mass_optimiz(p_mid2,q2,Me(2));
        
        L1 = length_ratio*R0*lhat1(i);
        L2 = length_ratio*R0*lhat2(i);
        
        p_mid2 = p_terminal + 8*mu*(L2/2)*q2/(pi*(lhat2(i)*R0)^4);
        p_bifur = p_terminal + 8*mu*(L2)*q2/(pi*(lhat2(i)*R0)^4);
        p_mid1 = p_bifur - 8*mu*(L1/2)*q1/(pi*(lhat1(i)*R0)^4);
        
        Y=[p_mid1,p_bifur,p_mid2];
        
        err = norm(abs(Y-Y0));
        
        Y0 = Y;
        c = c+1;
        
    end
    
end

y1 = log(Q_parent);
x1 = log(lhat1*R0);
P1 = polyfit(x1,y1,1);

P1(1);

y2 = log(Q_parent/2);
x2 = log(lhat2*R0);
P2 = polyfit(x2,y2,1);

P2(1);

