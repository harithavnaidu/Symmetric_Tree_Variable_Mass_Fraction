% modified by Vasilina, August 2, 2018
function [Mt,Me,SMCtoCOL,R,L,Table,p_mid] = TreeOptimization2(q_parent,p_terminal,N_gen) 
% This function solves the optimization problem of the arterial tree
% Given: q_parent,p_terminal,N_gen
% Output:

%% initialization
global R0 mu 

% preallocate arrays
Mt = zeros(1,N_gen);
Me = zeros(1,N_gen);
Mtotal = zeros(1,N_gen);
SMCtoCOL = zeros(1,N_gen); 
phi_e = zeros(1,N_gen); 

R = zeros(1,N_gen);
L = zeros(1,N_gen);
p_mid = zeros(1,N_gen); 
Res = zeros(1,N_gen);   
Table = zeros(N_gen,4);

% initialization
lhat = ones(1,N_gen);

% intitalize R,L,HydRes and Table
for k=1:N_gen
    R(k) = R0*lhat(k);
    L(k) = LengthSegmentk(R(k));
    Res(k) = 8*mu*(L(k))/(pi*R(k)^4);
end

% %initialization of hemodynamics
% Table = SymmetricArterialTree2(N_gen,q_parent,p_terminal,Res); 
[q, p_term]= SymmetricArterialTree2(N_gen,q_parent,p_terminal,Res); 

% p_term = Table(:,3);
Y0 = p_term;

%% Solve for terminal pressure to fit mass and geometry relations
% for lhat,R,L,Table

err = 1000;
c = 1;
while err > 1e-6
    for j=1:N_gen
        % initialization for each segment (generation)
        q_seg = q(j); %Table(j,2);
%         Res(j)= Table(j,4);
%         p_term(j) = Table(j,3);
        
        % p_term(j)=Table(j,3);
        % p_inp(j)= p_term(j)+Res(j)*q_seg
        % p_mid(j)= (p_term(j)+p_inp(j))=2*p_term(j)/2+Res(j)*q_seg/2
        %average pressure in a segment
        p_mid(j) = p_term(j)+ Res(j)*q_seg/2; 

        % optimization of the mass at each segment
        [Mtotal(j), lhat(j), phi_e(j), SMCtoCOL(j)] = mass_optimiz_tree(p_mid(j),q_seg,j);   
        Mt(j) = (1 - phi_e(j))*Mtotal(j);
        Me(j) = phi_e(j)*Mtotal(j);

        % update geometry R,L
        R(j) = R0*lhat(j);
        L(j) = LengthSegmentk(R(j));

        % update resistance
        Res(j) = 8*mu*(L(j))/(pi* R(j)^4); 
%         Table(j,4) = Res(j);

    end
    % update steady hemodynamics for entire tree
    % with new resitance
    [q, p_term] = SymmetricArterialTree2(N_gen,q_parent,p_terminal,Res); 

    % update terminal pressure
    Y = p_term;

    % compute residual error
    err = norm(abs(Y-Y0));

    % update reference solution
    Y0 = Y;
    c = c+1;
    
    % for output
    for j=1:N_gen
        Table(j,:) = [j, q(j), p_term(j), Res(j)];
    end
end

disp(['Tree Pressure Optimization iterations # ',num2str(c),', tollerance 1e-6']);
