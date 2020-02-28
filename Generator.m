function [towInteraction,MGen] = Generator(secondOrder,GInput,MInput)

%% ========================================================================
% Pre-Calculations / Pre-allocation
res = GInput.l/GInput.gridSpacing;
towVol = GInput.l*GInput.w*GInput.Tt;
RVEVol = GInput.L*GInput.Wp*GInput.RVEThickness;
Nt = RVEVol/towVol;
Nlayer = floor(Nt*GInput.l*GInput.w/(GInput.L^2));
Nt = round(RVEVol/towVol);
Nt1 = Nt;
MaxRealLength = Nt1*GInput.l;
RatioExtra = (GInput.L*GInput.Wp)/((GInput.L+GInput.l/4)*(GInput.Wp+GInput.l/4));
Nt=Nt/RatioExtra;
x=0:GInput.gridSpacing:GInput.Wp;
y=0:GInput.gridSpacing:GInput.L;
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
    [X00(i),Y00(i),angle(i),xv(i,1:5),yv(i,1:5),center1{i},center2{i},center3{i},center4{i},center5{i},topEdge{i},botEdge{i},Lef]=tow_placement(GInput.L,GInput.Wp,GInput.l,GInput.w,optAngles(i),res,GInput.gridSpacing);
    
    condition = floor(GInput.Vf_thresh/MInput.VfL*mean(t_counter))+2;
 
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

%% Thickness distribution of each tow & Volume fraction map!
for i=1:length(y)
    for j=1:length(x)
        aux2(:)=Cube3D(i,j,:);
        towPerNode{i,j}=aux2(aux2~=0);
        towCountPerNode(i,j)=length(aux2(aux2~=0));
        nodalPlyThickness(i,j) = 1/towCountPerNode(i,j)*GInput.RVEThickness;
        VFmap(i,j) = MInput.VfL * GInput.Tt/nodalPlyThickness(i,j);
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
    aux = interp2(temp,Center(:,1)/GInput.gridSpacing,Center(:,2)/GInput.gridSpacing);
    auxx=aux;
    auxx(isnan(auxx))=[];
    %numerical approximation
    if length(auxx)<=length(aux)/5
        aux(isnan(aux))=0;
    end
    waviness1{i} = aux;
    Center=center2{i};
    aux = interp2(temp,Center(:,1)/GInput.gridSpacing,Center(:,2)/GInput.gridSpacing);
    auxx=aux;
    auxx(isnan(auxx))=[];
    %numerical approximation
    if length(auxx)<=length(aux)/5
        aux(isnan(aux))=0;
    end
    waviness2{i} = aux;
    Center=center3{i};
    aux = interp2(temp,Center(:,1)/GInput.gridSpacing,Center(:,2)/GInput.gridSpacing);
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
        if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=GInput.L))~=0
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
            if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=GInput.L))~=0
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
            if length(aux10(aux10 <= -0.5))~=0 || length(aux10(aux10>=GInput.L))~=0
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
        if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=GInput.Wp))~=0
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
            if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=GInput.Wp))~=0
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
            if length(auxW(auxW <= -0.5))~=0 || length(auxW(auxW>=GInput.Wp))~=0
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
    xcoord = round(auxW1{i}(:,1)/GInput.gridSpacing);
    xcoord(xcoord<=0)=1;
    ycoord = round(aux9{i}(:,1)/GInput.gridSpacing);
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
    if max(aux7)>GInput.RVEThickness
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
    towInteraction.(temp).avgVf = MInput.VfT * GInput.Tt/towInteraction.(temp).avgThick;
end

MGen.res = res;
MGen.Nt = Nt;
MGen.Leffective = Leffective;
MGen.centerThickness = centerThickness;
MGen.RealWaviness = RealWaviness;
MGen.angle = angle;
MGen.VFmap = VFmap;



end

