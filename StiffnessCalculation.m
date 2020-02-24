function [Ew,MGen] = StiffnessCalculation(towInteraction,GInput,MInput,MGen)

%% ========================================================================
% ================== Tow Waviness Stiffness effect ========================
% =========================================================================

resol=MGen.res-2;
for i=1:MGen.Nt
    towLen{i} = linspace(0,MGen.Leffective(i),length(MGen.RealWaviness{i}));
end


% =========================================================================
Cwaviness = zeros(6);
for i = 1:floor(MGen.Nt)
    temp = strcat('towID_',num2str(i));
    VfTi = towInteraction.(temp).avgVf;
    VfL=VfTi;
    avgTt = mean(MGen.centerThickness{i});
    [C] = mesoSL_HT( MInput.Em, MInput.Gm, MInput.num,MInput.phi_f, GInput.l ,avgTt,GInput.w, MInput.VfL,VfTi, MInput.E11t,MInput.E22t,MInput.G12t,MInput.G23t,MInput.nu12t,MInput.nu23t );
    S = inv(C);
    % =====================================================================
    [Swaved,towOPAngles{i},ELongitudinal{i},Filtered{i}] = stiffnessAvgFFT(S,towLen{i},MGen.RealWaviness{i},GInput.l);
    towInteraction.(temp).outPlaneangles = towOPAngles{i};
    Savg = inPlaneRot(Swaved,MGen.angle(i));
    %
    temp = strcat('towID_',num2str(i));
    towInteraction.(temp).towAvgStiff = 1/Savg(1,1);
    towInteraction.(temp).towSTensorGlobal = Swaved;
    Cavg = inv(Savg);
    Cwaviness = Cwaviness + Cavg*MGen.Leffective(i);
end
Cwaviness = Cwaviness/sum(MGen.Leffective);
MGen.towOPAngles = towOPAngles;

Cnonwaviness=zeros(6);
for i = 1:floor(MGen.Nt)
    temp = strcat('towID_',num2str(i));
    VfTi = towInteraction.(temp).avgVf;
    VfL=VfTi;
    avgTt = mean(MGen.centerThickness{i});
    [C] = mesoSL_HT( MInput.Em, MInput.Gm, MInput.num,MInput.phi_f, GInput.l ,avgTt,GInput.w, MInput.VfL,VfTi, MInput.E11t,MInput.E22t,MInput.G12t,MInput.G23t,MInput.nu12t,MInput.nu23t );
    S = inv(C);
    Saux =inPlaneRot(S,MGen.angle(i));
    Caux = inv(Saux);
    Cnonwaviness = Cnonwaviness + Caux*MGen.Leffective(i);
end
Cnonwaviness = Cnonwaviness/sum(MGen.Leffective);

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


fprintf('E1 without waviness effects is %d.\n',Enw(1));
fprintf('E2 without waviness effects is %d.\n',Enw(2));
fprintf('============================================\n');
fprintf('E1 with waviness effects is %d.\n',Ew(1));
fprintf('E2 with waviness effects is %d.\n',Ew(2));
fprintf('============================================\n');
fprintf('E1 difference in percentage is %d.\n',percent1);
fprintf('E2 difference in percentage is %d.\n',percent2);
end

