% Halpin-Tsai Approach (Discontinuous Fibres) (Adopted from 'Analysis and Performance of Fibre Composites - Pg 141 by Bhagwan D.Agarwal' & 'The Halpin-Tsai Equations: A Review - 1976')
% Note: Halpin-Tsai does not differentiate between ellipsoidal and rectangular bundles
function [C] = mesoSL_HT( Em, Gm, nu_m, phi_f, l,t,w, Vfl,VfT, Et11,Et22,Gt12,Gt23,nu_t12,nu_t23 )
%% geometry
tchar=w*t/(2*(w+t));

% tchar=t/2;
lchar=l/8; 

% phi_t=sqrt(4*w*t/pi);
tm=(sqrt(pi/(4*Vfl))-1)*phi_f;

%% SLM for longitudinal stiffness
lambda=sqrt(2*Gm/(tchar*tm*Et11));
Et11 = Vfl/0.58*Et11;
El11=Et11/(1+1/(lambda*lchar*tanh(lambda*lchar)));
% El11 = Et11;

%% Halpin-Tsai Constants: Zeta and Eta %%
sl_t = sqrt(w*t);
% Zeta_L = 2*(l/sl_t);
% Eta_L = (Et11/Em - 1)/(Et11/Em + Zeta_L);
% El11 = Em*(1 + Zeta_L*Eta_L*Vtl)/(1 - Eta_L*Vtl); 

%%           Ellipsoidal/Rectangular/Circular/Square Chip

Zeta_T = 2*(w/t);
Zeta_D = 2*(t/w);                                                % Zeta_D for Out-of-Plane Modulus (E33 Hypothesised Extension)
Zeta_G12 = 1;                                                    % Zeta_G12 = 1 (Conservative)
Zeta_G23 = 1; 

Eta_T = (Et22/Em - 1)/(Et22/Em + Zeta_T);
if (w == t)
    Eta_D = Eta_T;                                               % Eta_D for Out-of-Plane Modulus (E33 Hypothesised Extension)
else
    Eta_D = (Et22/Em - 1)/(Et22/Em + Zeta_D);                    % Ellipsoidal/Rectangular Bundle: Eta_D for Out-of-Plane Modulus (Hypothesised Extension)
end

Eta_G12 = (Gt12/Gm - 1)/(Gt12/Gm + Zeta_G12);
Eta_G23 = (Gt23/Gm - 1)/(Gt23/Gm + Zeta_G23);

%% Lamina Properties %%
El22 = Em*(1 + Zeta_T*Eta_T*VfT)/(1 - Eta_T*VfT);                % Lamina Transverse Modulus
El33 = Em*(1 + Zeta_D*Eta_D*VfT)/(1 - Eta_D*VfT);                % Lamina Out-of-Plane Modulus (E33 Hypothesised Extension)

Gl12 = Gm*(1 + Zeta_G12*Eta_G12*VfT)/(1 - Eta_G12*VfT);          % Lamina In-Plane Shear Modulus
Gl23 = Gm*(1 + Zeta_G23*Eta_G23*VfT)/(1 - Eta_G23*VfT);          % Lamina Through-Thickness Shear Modulus
Gl13 = Gl12;

nu_l12 = nu_t12*VfT + nu_m*(1 - VfT);                            % Lamina Major Poisson Ratio
nu_l21 = nu_l12*El22/El11;                                       % Lamina Minor Poisson Ratio (Hooke's Law: Transversely Isotropic Bundle)

nu_l23 = nu_t23*VfT + nu_m*(1 - VfT);                            % Lamina Major Poisson Ratio
nu_l32 = nu_l23*El33/El22;     

nu_l13 = nu_l12;                                                 % Lamina Major Poisson Ratio
nu_l31 = nu_l21;     

%% Classic Laminate Theory Inputs %%
% S=[1/El11, -nu_l21/El22, -nu_l31/El33, 0, 0, 0;
%     -nu_l12/El11, 1/El22, -nu_l32/El33, 0, 0, 0;
%     -nu_l13/El11, -nu_l23/El22, 1/El33, 0, 0, 0;
%     0, 0, 0, 1/Gl23, 0, 0;
%     0, 0, 0, 0, 1/Gl13, 0;
%     0, 0, 0, 0, 0, 1/Gl12];

S=[1/El11, -nu_l21/El22, -nu_l31/El33, 0, 0, 0;
    -nu_l21/El22, 1/El22, -nu_l32/El33, 0, 0, 0;
    -nu_l31/El33, -nu_l32/El33, 1/El33, 0, 0, 0;
    0, 0, 0, 1/Gl23, 0, 0;
    0, 0, 0, 0, 1/Gl13, 0;
    0, 0, 0, 0, 0, 1/Gl12];
C=inv(S);

end