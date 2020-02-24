function [] = Plots(VF,L,Wp,towOPAngles)

oS = size(VF) ;
VFmap = VF(5:end-5,5:end-5);
VFmap = imresize(VFmap,oS);
%============== Volume fraction Variation and Distribution ================
%%
figure()
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'tex');
imagesc(VFmap);
% h = colorbar('horizontal');
h = colorbar;
ylabel(h, '$V^{\mathrm{f}}$')
h.Label.Interpreter = 'latex';
colormap jet
axis([0,Wp,0,L])
xlabel ('$x$(mm)')
ylabel ('$y$(mm)')
xticks([1 L+1])
xticklabels({0','$L$',})
yticks([1 Wp+1])
yticklabels({0','$W$',})

%% ========================================================================
% ==================== Volume fraction Distribution =======================
% =========================================================================
myFit = fitdist(VFmap(:), 'kernel');
index = linspace(min(VFmap(:)), max(VFmap(:)), 1000);
%
figure()
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'tex');
plot(index, pdf(myFit, index))
xlabel '$V^{\mathrm{f}}$'
ylabel 'Pr'
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',9,...
    'FontName','Times',...
    'XTick',[0:0.2:1],...
    'XLim',[0 1],...
    'YLim',[0 1.1*max(pdf(myFit, index))],...
    'Box','off');
box off

%% ========================================================================
% ======================= Tow angles Distribution =========================
% =========================================================================
towAnglesLin =cell2mat(towOPAngles);
myFit = fitdist(towAnglesLin'*180/pi, 'kernel');
index = linspace(min(towAnglesLin*180/pi), max(towAnglesLin*180/pi), 1000);

%
figure()
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter', 'tex');
plot(index, pdf(myFit, index))
xlabel '$\phi~^\circ$'
ylabel 'Pr'
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',9,...
    'FontName','Times',...
    'XTick',[-40:10:40 ],...
    'XLim',[-40 40],...
    'YLim',[0 1.1*max(pdf(myFit, index))],...
    'Box','off');
box off

end

