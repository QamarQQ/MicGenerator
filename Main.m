% Copyright 2018 Marco Alves

% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at

%     http://www.apache.org/licenses/LICENSE-2.0

% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%% ========================================================================
% Authors: Marco Alves    email: marco.alves@imperial.ac.uk
%          Soraia Pimenta emai: soraia.pimenta@imperial.ac.uk
% Imperial College London - meComposites

%%
clear
clc
close all
%% 
% Define materialinputs
MInput.VfL = 0.57;      % fibre volume fraction in composites 
MInput.VtL = 1;          % tow volume fraction in composites
MInput.VfT = MInput.VfL/MInput.VtL;    % fibre volume fraction in tow
MInput.Vmt = 1 - MInput.VfT;    % matrix volume fraction in tow
MInput.E11t=116e3;
MInput.E22t=9e3;
MInput.Em=4e3;
MInput.G12t=5.6e3;
MInput.G23t=2800;
MInput.Gm=1.0e3;
MInput.nu12t=0.34;
MInput.nu23t=0.34;
MInput.num=0.4;
MInput.phi_f=0.007; % fibre diameter

% Define the geometry inputs;
GInput.l = 50;                           % tow length
GInput.w = 8;                            % tow width
GInput.Tt=0.1;                           % tow thickness 
GInput.L = 150;                          % length of the specimen
GInput.Wp = 150;                         % width of the specimen
GInput.RVEThickness = 3;                 % thickness of the specimen
GInput.gridSpacing=1;                    % grid Spacing value
GInput.Vf_thresh = 0.75;                 % volume fraction threshold value

% Planar second order tensor describing the overall tow orientation state
secondOrder=[0.5  , 0  ;0   ,0.5];

%%
[towInteraction,MGen] = Generator(secondOrder,GInput,MInput);

[Ew,MGen] = StiffnessCalculation(towInteraction,GInput,MInput,MGen);

Plots(MGen.VFmap,GInput.L,GInput.Wp,MGen.towOPAngles)


