function [ Savg,towAngle,E1Longitudinal,Yfiltered ] = stiffnessAvgFFT( S,towLen,realWaviness,towLength,gridSpacing )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Data = realWaviness;
N = towLen;

L=towLen(end);    
res=length(Data);
dx = L/res;               %spatial discretization
m = mean(Data);
Data = Data-m;
n = length(Data);
N = L*linspace(0,1,n);
Ks = 1/dx;                %spatial sampling frequency
k = Ks*(-n/2:n/2-1)/n;

% Window function - low pass filter: kc = cut off freq. -> 1/kc = lambda (filtered)  -|-  Lk = length scale of tanh -> steepness of the transition
%kc = 0.1176; % This is filtering lengths = tow length
% kc=1/gridSpacing;
kc=0.3;
%
Lk = 0.1;
window = (tanh((k+kc)./Lk) - tanh((k-kc)./Lk))/2;
%
Y = fft(Data);
Yfiltered = real(ifft(ifftshift(window.*fftshift(Y))));
%
% figure()
% plot(realWaviness)
% hold on
% plot(Yfiltered)
%
towAngle = atan(diff(Yfiltered)./diff(N));
towAngleDegree = towAngle*180/pi;

for i=1:length(towLen)-1
    incrementLength(i) = (towLen(i+1)-towLen(i))/(cos(towAngle(i)));
end


Cavg=zeros(6);
for i = 1:length(towLen)-1
    Saux = outPlaneRot(S,towAngle(i));
    Caux = inv(Saux);
    E1Longitudinal(i) = 1/Saux(1,1);
    Cavg = Cavg + Caux*incrementLength(i);
end

Cavg = Cavg/(sum(incrementLength));
Savg = inv(Cavg);



end



