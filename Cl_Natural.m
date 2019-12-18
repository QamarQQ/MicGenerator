function [ a4 ] = Natural( a2 )
%Natural Natural closure approximation
%   Natural closure approximation to obtain the 4th order orientation
%   tensor from the 2nd order tensor
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    a4(i,j,k,l)=1/6*det(a2)*(kroneckerDelta(i,j)*kroneckerDelta(k,l)+kroneckerDelta(i,k)*kroneckerDelta(j,l)+kroneckerDelta(i,l)*kroneckerDelta(j,k))+1/3*(a2(i,j)*a2(k,l)+a2(i,k)*a2(j,l)+a2(i,l)+a2(j,k));
                end
            end
        end
    end
end

