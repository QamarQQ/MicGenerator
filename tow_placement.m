%FUNCTION. Generate tow centre's xy positions and orientation randomly

function [X00,Y00,angle,xv,yv,center1,center2,center3,center4,center5,topEdge,botEdge,Lef]=tow_placement(L,Wp,l,w,angle,res,gridSpacing)

%     X00=gridSpacing + (Wp-gridSpacing)*rand();   %x,y coord of tow's centre pt lies within 0 to RVE_size
%     Y00=gridSpacing + (L-gridSpacing)*rand();
    X00= (Wp)*rand();   %x,y coord of tow's centre pt lies within 0 to RVE_size
    Y00= (L)*rand();
    %tow orientation
%     angle=rand()*pi;        %from 0to pi, -> 0to180deg
    x2 = [-w/2  w/2 w/2 -w/2 -w/2];     %x2, y2 used to compute corner coordinates of a tow, based on tow width and length 
    y2 = [-l/2 -l/2 l/2 l/2 -l/2];  
    %record corner coordinates of a tow. Got 5 pts (from bottom left corner, goes c-clockwise to bottom left again)
    xv=[X00 + cos(angle)*x2 - sin(angle)*y2]';  %5 pts outlining the x-coord of corners of rectangle, from bottom left, goes c-clockwise to bottom left 
    yv=[Y00 + sin(angle)*x2 + cos(angle)*y2]';  %5 pts outlining the y-coord of corners of rectangle, same:c-clockwise
    
centerPx(1)=(xv(1)+xv(2))/2;
centerPy(1)=(yv(1)+yv(2))/2;
centerPx(2)=(xv(3)+xv(4))/2;
centerPy(2)=(yv(3)+yv(4))/2;
center=[linspace(centerPx(1),centerPx(2),res).' linspace(centerPy(1),centerPy(2),res).'];
Leff = center;
[row,col]=size(center);

% Save the edge locations of the tows
topedge(:,1) = linspace(xv(1),xv(2),l/(gridSpacing/2));
topedge(:,2) = linspace(yv(1),yv(2),l/(gridSpacing/2)) ;
botedge(:,1) = linspace(xv(3),xv(4),l/(gridSpacing/2)) ;
botedge(:,2) = linspace(yv(3),yv(4),l/(gridSpacing/2)) ;

topEdge(:,1) = (round(topedge(:,1)./gridSpacing));
topEdge(:,2) = (round(topedge(:,2)./gridSpacing));
botEdge(:,1) = (round(botedge(:,1)./gridSpacing));
botEdge(:,2) = (round(botedge(:,2)./gridSpacing));

% offSet = 0.1*res;
offSet = 0.00*res;
center1=center(1+offSet:row-offSet,:);

midleftx(1)=(xv(2)+centerPx(1))/2;
midlefty(1)=(yv(2)+centerPy(1))/2;
midleftx(2)=(xv(3)+centerPx(2))/2;
midlefty(2)=(yv(3)+centerPy(2))/2;
midl=[linspace(midleftx(1),midleftx(2),res).' linspace(midlefty(1),midlefty(2),res).'];
center2=midl(1+offSet:row-offSet,:);

midrightx(1)=(xv(1)+centerPx(1))/2;
midrighty(1)=(yv(1)+centerPy(1))/2;
midrightx(2)=(xv(4)+centerPx(2))/2;
midrighty(2)=(yv(4)+centerPy(2))/2;
midr=[linspace(midrightx(1),midrightx(2),res).' linspace(midrighty(1),midrighty(2),res).'];
center3=midr(1+offSet:row-offSet,:);


midl=[linspace(xv(1)+1,xv(4)+1,res).' linspace(yv(1),yv(4),res).'];
center4=midl(1+offSet:row-offSet,:);

midr=[linspace(xv(2)-1,xv(3)-1,res).' linspace(yv(2),yv(3),res).'];
center5=midr(1+offSet:row-offSet,:);


indices1 = find(Leff(:,1)<0);
indices2 = find(Leff(:,2)<0);
%
indices3 = find(Leff(:,1)>Wp);
indices4 = find(Leff(:,2)>L);

catIndex = [indices1; indices2; indices3; indices4];
catIndex = unique(sort(catIndex));
%
Leff(catIndex,:) = [];

if size(Leff,1) == 1
    Lef = gridSpacing;
else
    Lef = pdist([Leff(1,1) Leff(1,2); Leff(end,1) Leff(end,2)],'euclidean');
end
end