function[E11t, E22t, G12t, G23t, nu12t, nu23t, Em,Gm, num,...
    VfL, VtL, VfT, Vmt,phi_f] = getInputs()

VfL = 0.57;      % fibre volume fraction in composites 
VtL = 1;          % tow volume fraction in composites
VfT = VfL/VtL;    % fibre volume fraction in tow
Vmt = 1 - VfT;    % matrix volume fraction in tow

%%
% E11t=124e3;
E11t=116e3;
E22t=9e3;
Em=4e3;
G12t=5.6e3;
G23t=2800;
Gm=1.0e3;
nu12t=0.34;
nu23t=0.34;
num=0.4;
% Lt=50;     % tow length
% wt=8;      % tow width
% Tt=0.164;    % tow thickness

phi_f=0.007; % fibre diameter


end