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
addpath('AuxiliaryFunctions')
%% INITIALIZATION
% Get the material inputs;
[E11t, E22t, G12t, G23t, nu12t, nu23t, Em,Gm, num, VfL, VtL, VfT, Vmt, phi_f] = getInputs();

% Define the geometry inputs;
l = 50;                           % tow length
w = 8;                            % tow width
Tt=0.1;                           % tow thickness 
L = 150;                          % length of the specimen
Wp = 150;                         % width of the specimen
RVEThickness = 1;                 % thickness of the specimen
gridSpacing=1;                    % grid Spacing value
Vf_thresh = 0.75;                 % volume fraction threshold value

% Planar second order tensor describing the overall tow orientation state
secondOrder=[0.5  , 0  ;0   ,0.5];
%% ========================================================================
% Pre-Calculations / Pre-allocation
res = l/gridSpacing;
TolArea = 0.05;
towVol = l*w*Tt; 
RVEVol = L*Wp*RVEThickness;
Nt = RVEVol/towVol;
Nlayer = floor(Nt*l*w/(L^2));
Nt = round(RVEVol/towVol);
Nt1 = Nt;
MaxRealLength = Nt1*l;
RatioExtra = (L*Wp)/((L+l/4)*(Wp+l/4));
Nt=Nt/RatioExtra;
x=0:gridSpacing:Wp;
y=0:gridSpacing:L;
Ngridx = length(x);
Ngridy = length(y);
[X,Y] = meshgrid(x,y);
t_counter_x=X(:)';
t_counter_y=Y(:)';
mCount = zeros(length(y),length(x));
mCount1 = mCount(:);

%% GENERATE ANGLES ACOORDING TO SECOND ORDER ORIENTATION TENSOR

[psi,angle]=FOD(secondOrder);
aux1=2*psi*Nt*pi/2;
sum(aux1);
aux1=round(aux1);
optAngles(1:aux1(1))=pi/180;
temp=length(optAngles);
for i=2:length(aux1)
    optAngles(temp+1:temp+aux1(i))=i*pi/180;  
    temp=length(optAngles);
end
dif=Nt-length(optAngles);
aux=length(optAngles);
for i=1:dif
    optAngles(aux+i)=optAngles(randi(aux));
end
optAngles=optAngles(randperm(length(optAngles)));

%% PLACE TOWS IN RVE

figure()
set(0,'DefaultTextInterpreter', 'latex')
xlabel ('$x$(mm)')
ylabel ('$y$(mm)')
axis equal
axis ([0 250 0 250])
set(gca,'TickLabelInterpreter','latex')
hold on
plot(t_counter_x,t_counter_y,'r+','MarkerSize',5)
hold on
t_counter = zeros(1,length(t_counter_x));
t_counter0 = t_counter;

Cube3D=zeros(length(y),length(x),20*Nlayer);

EffectiveLengthTot = 0;
i=1;

while EffectiveLengthTot<=MaxRealLength
    %Random tow placement
    [X00(i),Y00(i),angle(i),xv(i,1:5),yv(i,1:5),center1{i},center2{i},center3{i},center4{i},center5{i},topEdge{i},botEdge{i},Lef]=tow_placement(L,Wp,l,w,optAngles(i),res,gridSpacing);
    
    condition = floor(Vf_thresh/VfL*mean(t_counter))+2;
 
    t_counter = t_counter + inpolygon(t_counter_x, t_counter_y, xv(i,:), yv(i,:));
    [in,on]=inpolygon(t_counter_x, t_counter_y, xv(i,:), yv(i,:));
%     if max(t_counter)-min(t_counter)> Delta_l(i)
    if max(t_counter)>condition
        t_counter = t_counter0; 
    else
        
        patch(xv(i,:),yv(i,:),'b');
        mCount1(in)=mCount1(in)+1;
        index=find(in~=0);
        towPos{i}= in;
        towCount{i}= (reshape(mCount1,Ngridy,Ngridx));
        towPosition{i}=mCount1(in);
        positioningTow=i;
        t_counter0 = t_counter;
        % ============= Waviness information =============
        aux1=zeros(length(t_counter_y),1);
        aux=towCount{i};
%         aux=flip(aux);
        in1=towPos{i};
        aux=aux(:);
        aux=aux(in1);
        aux1(in1)=aux;
        towWave{i}=(reshape(aux1,Ngridy,Ngridx));
        % ============== Update the 3D cube ==============
        aux1=flip(towWave{i});
        maxIndex=max(max(aux1));
        minIndex=min(min(aux1(aux1>0)));
        nLayers=[minIndex:maxIndex];
        for k = 1:length(nLayers)
            layerIndex{k} = find(aux1==nLayers(k));
            pos=zeros(length(y),length(x));
            pos(layerIndex{k})=i;
            pos=flip(pos);
            Cube3D(:,:,nLayers(k))=Cube3D(:,:,nLayers(k)) + pos;
        end
        i=i+1;
        EffectiveLengthTot = EffectiveLengthTot + Lef;
    end
    i
end
Nt = i-1;

close all

%% ========================================================================
% ================= Individual tow interaction build ======================
% =========================================================================

for i=1:Nt
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp) = TowInteractions;
end

for i = 1:Nt
    % Assign center line to each tow
    top = zeros(length(x),length(y));
    bot = zeros(length(x),length(y));
    [r,c,v] = ind2sub(size(Cube3D),find(Cube3D == i));
    for k = 1:length(r);
        top(r(k),c(k)) = Cube3D(r(k),c(k),v(k)+1);
        if v(k) <= 1
            bot(r(k),c(k)) = 0;
        else
            bot(r(k),c(k)) = Cube3D(r(k),c(k),v(k)-1);
        end
    end
    top1 = flip(top);
    bot1 = flip(bot);
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp).topInteractionsMap = top1;
    towInteraction.(temp).bottomInteractionsMap = bot1;
    
    %   ======================== bottom interactions ======================
    [b,m1,n1] = unique(bot,'first');
    [c1,d1] =sort(m1);
    botInteract = b(d1);
    botInteract(botInteract==0)=[];
    if isempty(botInteract) == 1
        botInteractArea=[];
    else
        for j=1:length(botInteract)
            botInteractArea(j) = length(bot(bot==botInteract(j)))/length(bot(bot>0));
        end
    end
    
    invalidBottomInteractions = botInteract(find(botInteractArea <= TolArea));
    botInteract = botInteract(find(botInteractArea > TolArea));
    botInteractArea = botInteractArea(find(botInteractArea>TolArea));
    
    % Remove invalid interactions from the interaction map ===
    for j = 1:length(invalidBottomInteractions)
        bot(find(bot==invalidBottomInteractions(j))) = 0;
    end
    bot1 = flip(bot);
    towInteraction.(temp).bottomValidInteractionsMap = bot1;
    towInteraction.(temp).bottomInteractingTows = botInteract';
    % 
    towInteraction.(temp).bottomInteractionAngles = optAngles(botInteract);
    towInteraction.(temp).bottomInteractionAreas = botInteractArea;
    
    %    ======================== top  interactions =======================
    [b,m1,n1] = unique(top,'first');
    [c1,d1] =sort(m1);
    topInteract = b(d1);
    topInteract(topInteract==0)=[];
    for j=1:length(topInteract)
        topInteractArea(j) = length(top(top==topInteract(j)))/length(top(top>0));
    end
    
    if isempty(topInteract) == 1
        1;
    else
        invalidTopInteractions = topInteract(find(topInteractArea <= TolArea));
        topInteract = topInteract(find(topInteractArea > TolArea));
        topInteractArea = topInteractArea(find(topInteractArea>TolArea));
    end
    % Remove invalid interactions from the interaction map ===
    for j = 1:length(invalidTopInteractions)
        top(find(top==invalidTopInteractions(j))) = 0;
    end
    top1 = flip(top);
    towInteraction.(temp).topValidInteractionsMap = top1;
    towInteraction.(temp).topInteractingTows = topInteract';
    %    
    towInteraction.(temp).topInteractionAngles = optAngles(topInteract);
    towInteraction.(temp).topInteractionAreas = topInteractArea;
    towInteraction.(temp).totalInteractions = length(topInteract) + length(botInteract);
    towInteraction.(temp).towAngle = optAngles(i)';
    % Reset variables ======
    topInteract = [];
    topInteractArea = [];
    botInteract = [];
    botInteractArea = [];
end
% % =========================================================================
% % =========================================================================
% % Thickness distribution of each tow & Volume fraction map!
for i=1:length(y)
    for j=1:length(x)
        aux2(:)=Cube3D(i,j,:);
        towPerNode{i,j}=aux2(aux2~=0);
        towCountPerNode(i,j)=length(aux2(aux2~=0));
        nodalPlyThickness(i,j) = 1/towCountPerNode(i,j)*RVEThickness;
        VFmap(i,j) = VfL * Tt/nodalPlyThickness(i,j);
    end
end
nodalPlyThickness(nodalPlyThickness==Inf) = 0;

%% ========================================================================
% ===================== Tow Waviness Description ==========================
% =========================================================================

% Layers that each tow is positioned into
for i =1:floor(Nt)
    aux1=zeros(length(t_counter_y),1);
    aux=towCount{i};
    in=towPos{i};
    aux=aux(:);
    aux=aux(in);
    aux1(in)=aux;
    towWave{i}=(reshape(aux1,Ngridy,Ngridx));
end

% Interpolating the waviness along the axial tow direction
for i =1:floor(Nt)
    Center=center1{i};
    temp=towWave{i};
    aux = interp2(temp,Center(:,1)/gridSpacing,Center(:,2)/gridSpacing);
    auxx=aux;
    auxx(isnan(auxx))=[];
    %numerical approximation
    if length(auxx)<=length(aux)/5
        aux(isnan(aux))=0;
    end
    waviness1{i} = aux;
    Center=center2{i};
    aux = interp2(temp,Center(:,1)/gridSpacing,Center(:,2)/gridSpacing);
    auxx=aux;
    auxx(isnan(auxx))=[];
    %numerical approximation
    if length(auxx)<=length(aux)/5
        aux(isnan(aux))=0;
    end
    waviness2{i} = aux;
    Center=center3{i};
    aux = interp2(temp,Center(:,1)/gridSpacing,Center(:,2)/gridSpacing);
    auxx=aux;
    auxx(isnan(auxx))=[];
    %numerical approximation
    if length(auxx)<=length(aux)/5
        aux(isnan(aux))=0;
    end
    waviness3{i}=aux;
end

for i =1:floor(Nt)
    Center=center1{i};
    % Trimming the tows that fall out of the grid space ===================
    conta=0;
    for j =1:length(Center)
        aux10=Center(j,2);
        if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=L))~=0
            1;
        else
            conta=conta+1;
            aux9{i}(conta,:) = Center(j,2);
        end
    end
    if length(aux9{i}) < 4
        Center=center2{i};
        conta=0;
        for j =1:length(Center)
            aux10=Center(j,2);
            if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=L))~=0
                1
            else
                conta=conta+1;
                aux9{i}(conta,:) = Center(j,2);
            end
        end
    end
    if length(aux9{i})<4
        Center=center3{i};
        conta=0;
        for j =1:length(Center)
            aux10=Center(j,2);
            if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=L))~=0
                1;
            else
                conta=conta+1;
                aux9{i}(conta,:) = Center(j,2);
            end
        end
    end
    
    % Trim in x direction 
    
     conta1=0;
    for j =1:length(Center)
        auxW=Center(j,1);
        if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=Wp))~=0
            1;
        else
            conta1=conta1+1;
            auxW1{i}(conta1,:) = Center(j,1);
        end
    end
    if length(auxW1{i})<4
        Center=center2{i};
        conta1=0;
        for j =1:length(Center)
            auxW=Center(j,1);
            if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=Wp))~=0
                1;
            else
                conta1=conta1+1;
                auxW1{i}(conta1,:) = Center(j,1);
            end
        end
    end
    if length(auxW1{i})<4
        Center=center3{i};
        conta1=0;
        for j =1:length(Center)
            auxW=Center(j,1);
            if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=Wp))~=0
                1;
            else
                conta1=conta1+1;
                auxW1{i}(conta1,:) = Center(j,1);
            end
        end
    end
    temp = strcat('towID_',num2str(i));
    if length(aux9{i}) < length(auxW1{i})
        auxW1{i}(length(aux9{i})+1:end) = [];
    elseif length(aux9{i})  > length(auxW1{i})
        aux9{i}(length(auxW1{i})+1:end) = [];
    end
    towInteraction.(temp).centerLine = [auxW1{i} aux9{i}];
        
    
    
    %======================================================================
    Leffective(i) = pdist([auxW1{i}(1,1),aux9{i}(1); auxW1{i}(end),aux9{i}(end)],'euclidean');
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp).lEffective = Leffective(i);
    xcoord = round(auxW1{i}(:,1)/gridSpacing);
    xcoord(xcoord<=0)=1;
    ycoord = round(aux9{i}(:,1)/gridSpacing);
    ycoord(ycoord<=0)=1;
    for k=1:length(xcoord)
        centerThickness{i}(k)= nodalPlyThickness(ycoord(k),xcoord(k));
    end
    temp1=(waviness1{i}+waviness2{i}+waviness3{i})/3;
    temp1(isnan(temp1))=[];
    wavinessAvg{i}=temp1;
end

for i =1:floor(Nt)
    aux8 = wavinessAvg{i};
    aux7 = centerThickness{i};
    if max(aux7)>RVEThickness
        error('centerThickness cannot be bigger than RVEThickness')
    end
    if length(aux7)<length(aux8)
        aux8=aux8(1:length(aux7));
    elseif length(aux7)>length(aux8)
        aux7=aux7(1:length(aux8));
    end
    for k = 1:length(aux8)
        RealWaviness{i}(k) = aux7(k)*(aux8(k)-0.5);
    end
    %
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp).realThick = centerThickness{i};
    towInteraction.(temp).avgThick = mean(centerThickness{i});
    towInteraction.(temp).avgVf = VfT * Tt/towInteraction.(temp).avgThick;
end
%% ========================================================================
% ================== Tow Waviness Stiffness effect ========================
% =========================================================================

resol=res-2;
for i=1:Nt
    towLen{i} = linspace(0,Leffective(i),length(RealWaviness{i}));
end


% =========================================================================
Cwaviness = zeros(6);
for i = 1:floor(Nt)
    temp = strcat('towID_',num2str(i));
    VfTi = towInteraction.(temp).avgVf;
    VfL=VfTi;
    avgTt = mean(centerThickness{i});
    [C] = mesoSL_HT( Em, Gm, num,phi_f, l ,avgTt,w, VfL,VfTi, E11t,E22t,G12t,G23t,nu12t,nu23t );
    S = inv(C);
    % =====================================================================
    [Swaved,towOPAngles{i},ELongitudinal{i},Filtered{i}] = stiffnessAvgFFT(S,towLen{i},RealWaviness{i},l);
    towInteraction.(temp).outPlaneangles = towOPAngles{i};
    Savg = inPlaneRot(Swaved,angle(i));
    %
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp).towAvgStiff = 1/Savg(1,1);
    towInteraction.(temp).towSTensorGlobal = Swaved;
    Cavg = inv(Savg);
    Cwaviness = Cwaviness + Cavg*Leffective(i);
end
Cwaviness = Cwaviness/sum(Leffective);

Cnonwaviness=zeros(6);
for i = 1:floor(Nt)
    temp = strcat('towID_',num2str(i));
    VfTi = towInteraction.(temp).avgVf;
    VfL=VfTi;
    avgTt = mean(centerThickness{i});
    [C] = mesoSL_HT( Em, Gm, num,phi_f, l ,avgTt,w, VfL,VfTi, E11t,E22t,G12t,G23t,nu12t,nu23t );
    S = inv(C);
    Saux =inPlaneRot(S,angle(i));
    Caux = inv(Saux);
    Cnonwaviness = Cnonwaviness + Caux*Leffective(i);
end
Cnonwaviness = Cnonwaviness/sum(Leffective);

%=======  In plane modulus values  ======
Swaviness = inv(Cwaviness);
Snonwaviness = inv(Cnonwaviness);
%
Ew(1)=1/Swaviness(1,1);
Ew(2)=1/Swaviness(2,2);
%
Enw(1)=1/Snonwaviness(1,1);
Enw(2)=1/Snonwaviness(2,2);
%
percent1=(Enw(1)-Ew(1))/Enw(1)*100 ;
percent2=(Enw(2)-Ew(2))/Enw(2)*100 ;
%
fprintf('E1 difference in percentage is %d.\n',percent1);
fprintf('E2 difference in percentage is %d.\n',percent2);
fprintf('============================================\n');
fprintf('E1 without waviness effects is %d.\n',Enw(1));
fprintf('E2 without waviness effects is %d.\n',Enw(2));
fprintf('============================================\n');
fprintf('E1 with waviness effects is %d.\n',Ew(1));
fprintf('E2 with waviness effects is %d.\n',Ew(2));


%% ======================== Usefull Plots =========================

Plots(VFmap,L,Wp,towOPAngles)
