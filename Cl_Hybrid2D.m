function [ a4 ] = Hybrid2D( a2 )
    %Hybrid Hybrid closure approximation
    %   Hybrid closure approximation to obtain 4th order orientation tensor
    %   from 2nd order tensor
    % linear (planar) closure model
    a4_lin=Cl_Linear2D(a2);
    % quadratic closure
    a4_quad=Cl_Quadratic(a2);
    for i=1:2
        for j=1:2
            if i==1 && j==1
                f=2*a2(i,j)*a2(j,i)-1;
            else
                f=f+2*a2(i,j)*a2(j,i);
            end
        end
    end
    % f=1-27*det(a2);
    a4=(1-f)*a4_lin+f*a4_quad;
end

