%FUNCTION. Generate tow centre's xy positions and orientation randomly

function [xv,yv]=tow_placement(X00,Y00,L,l,w,angle)

    X00=L*rand();   %x,y coord of tow's centre pt lies within 0 to RVE_size
    Y00=L*rand();
    %tow orientation
%     angle=rand()*pi;        %from 0to pi, -> 0to180deg
    x2 = [-w/2  w/2 w/2 -w/2 -w/2];     %x2, y2 used to compute corner coordinates of a tow, based on tow width and length 
    y2 = [-l/2 -l/2 l/2 l/2 -l/2];  
    %record corner coordinates of a tow. Got 5 pts (from bottom left corner, goes c-clockwise to bottom left again)
    xv=[X00 + cos(angle)*x2 - sin(angle)*y2]';  %5 pts outlining the x-coord of corners of rectangle, from bottom left, goes c-clockwise to bottom left 
    yv=[Y00 + sin(angle)*x2 + cos(angle)*y2]';  %5 pts outlining the y-coord of corners of rectangle, same:c-clockwise
    
end