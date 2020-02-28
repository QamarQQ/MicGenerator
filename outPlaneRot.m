function [ Srot ] = outPlaneRot( S, angle )
% Marco Alves
% Out-of-plane rotation of tensor

m=cos(angle);
n=sin(angle);

T=[m^2   0   n^2   0   2*m*n  0;
    0    1    0    0     0    0;
   n^2   0   m^2   0  -2*m*n  0;
    0    0    0    m     0   -n;
  -m*n   0   m*n   0  m^2-n^2 0 ;
    0    0    0    n     0    m];

R=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 2 0 0;
   0 0 0 0 2 0;
   0 0 0 0 0 2];

Srot=R*inv(T)*inv(R)*S*T;

end

