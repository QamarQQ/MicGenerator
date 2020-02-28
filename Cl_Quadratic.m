function [ a4 ] = Quadratic( a2 )
%Quadratic Quadratic closure approximation
%   Quadratic closure approximation to find 4th order orientation tensor
%   from 2nd order tensor
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    a4(i,j,k,l)=a2(i,j)*a2(k,l);
                end
            end
        end
    end
end

