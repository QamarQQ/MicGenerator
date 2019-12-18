% Mori-Tanaka Approach %
% Inclusion sees the surrounding composite material by assuming matrix component in strain concentration tensor term

function [Cl] = mesoMT(Em, Gm, nu_m, l,t,w, Vtl, Et11,Et22,Gt12,nu_t12,nu_t23 )
    % Inputs for Chip Witdth/Thicknesss %%
     phi_t = sqrt(w*t*4/pi);
     
     %% tow stiffness tensor
     
     %% Compliance Tensor ('Fundamentals of Fibre Reinforced Composite Materials' A R Bunsell) %%
    S11 = 1/Et11;                          
    S22 = 1/Et22;
    S12 = - nu_t12/Et11;                                                            % S12 = -Yb12/Eb11 = -Yb21/Eb22
    S23 = - nu_t23/Et22;
    S66 = 1/Gt12;                                                                 

    
    % Bundle Compliance/Stiffness Tensor
    St = zeros(6);
    St(1,1) = S11;
    St(1,2) = S12;
    St(1,3) = S12;
    St(2,1) = S12;
    St(3,1) = S12;
    St(2,2) = S22;
    St(2,3) = S23;
    St(3,2) = S23;
    St(3,3) = S22;
    St(4,4) = (S22 - S23)/2;                                  
    St(5,5) = S66;
    St(6,6) = S66;
    Ct = St^(-1);
  

%% Matrix Stiffness Tensor ('Strength of Fibrous Composite' - Zheng Ming Huang) %%
    Sm = zeros(6);
    Sm(1,1) = 1/Em;
    Sm(1,2) = - nu_m/Em;
    Sm(1,3) = - nu_m/Em;
    Sm(2,1) = - nu_m/Em;
    Sm(2,2) = 1/Em;
    Sm(2,3) = - nu_m/Em;
    Sm(3,1) = - nu_m/Em;
    Sm(3,2) = - nu_m/Em;
    Sm(3,3) = 1/Em;
    Sm(4,4) = 1/Gm;
    Sm(5,5) = 1/Gm;
    Sm(6,6) = 1/Gm;
    Cm = Sm^(-1) ;


    %% Eshelby Terms for Prolate Spheroidal Shape in Isotropic Matrix (Adopted from Qiu and Weng 1990 Appendix A) %%            
    Alpha = l/phi_t;                                                                             % Aspect Ratio for Spherioidal Bundle/Inclusion
    g = (Alpha/(Alpha^2 - 1)^1.5)*(Alpha*(Alpha^2 - 1)^0.5 - acosh(Alpha));                                        % g for Prolate Shape for a > 1 (Different Expression for a<1 - 'Flow Induced Alignment in Composite Material Pg 302')
    S1111 = (1/(2 - 2*nu_m))*(1 - 2*nu_m + (3*Alpha^2 - 1)/(Alpha^2 - 1) - (1 - 2*nu_m + 3*Alpha^2/(Alpha^2 - 1))*g);
    S2222 = (3/(8 - 8*nu_m))*(Alpha^2/(Alpha^2 - 1)) + (1/(4 - 4*nu_m))*(1 - 2*nu_m - 9/(4*Alpha^2 - 4))*g;
    S2323 = (1/(4 - 4*nu_m))*(Alpha^2/(2*Alpha^2 - 2) + (1 - 2*nu_m - 3/(4*Alpha^2 - 4))*g);
    S1313 = (1/(4 - 4*nu_m))*(1 - 2*nu_m - (Alpha^2 + 1)/(Alpha^2 - 1) - 0.5*(1 - 2*nu_m - 3*(Alpha^2 + 1)/(Alpha^2 - 1))*g); 
    S1122 = (- 1/(2 - 2*nu_m))*(1 - 2*nu_m + 1/(Alpha^2 - 1)) + (1/(2 - 2*nu_m))*(1 - 2*nu_m + 3/(2*Alpha^2 - 2))*g;
    S2211 = (- 1/(2 - 2*nu_m))*(Alpha^2/(Alpha^2 - 1)) + (1/(4 - 4*nu_m))*((3*Alpha^2/(Alpha^2 - 1)) - (1 - 2*nu_m))*g;
    S2233 = (1/(4 - 4*nu_m))*((Alpha^2/(2*Alpha^2 - 2)) - (1 - 2*nu_m + 3/(4*Alpha^2 - 4))*g);

    Em_Eshelby = zeros(6);
    Em_Eshelby(1,1) = S1111;
    Em_Eshelby(2,2) = S2222;
    Em_Eshelby(3,3) = Em_Eshelby(2,2);
    Em_Eshelby(4,4) = 2*S2323;
    Em_Eshelby(5,5) = 2*S1313;
    Em_Eshelby(6,6) = Em_Eshelby(5,5);
    Em_Eshelby(1,2) = S1122;
    Em_Eshelby(1,3) = Em_Eshelby(1,2);
    Em_Eshelby(2,1) = S2211;
    Em_Eshelby(3,1) = Em_Eshelby(2,1);
    Em_Eshelby(2,3) = S2233;
    Em_Eshelby(3,2) = Em_Eshelby(2,3);


    %% Lamina Stiffness Tensor (Transversely Isotropic) based on Mori-Tanaka Equations - Tucker's Notation %%
    I = eye(6);                                                                                                           % 6x6 Unity Tensor
    A_Eshelby = (I + Em_Eshelby*Sm*(Ct - Cm))^(-1);                   % Eshelby Strain Concentrator
    A_MT = A_Eshelby*((1 - Vtl)*I + Vtl*A_Eshelby)^(-1);                                                                    % Mori-Tanaka Strain Concentrator
    Cl = Cm + Vtl*(Ct - Cm)*A_MT;                        % Lamina Stiffness Tensor


    %%  Classic Laminate Theory Inputs: C11, C22, C12, C66 %%
%     C11 = Cl(1,1) - (Cl(1,2)^2/Cl(2,2));
%     C22 = Cl(2,2) - (Cl(2,3)^2/Cl(2,2));
%     C12 = Cl(1,2) - Cl(1,2)*Cl(2,3)/Cl(2,2);
%     C66 = Cl(6,6);

    
end