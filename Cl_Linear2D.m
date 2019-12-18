function [ a4 ] = Linear2D( a2 )
%Linear Linear closure approximation
%   Linear closure approximation to obtain fourth order orientation tensor
%   from 2nd order tensor
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    a4(i,j,k,l)=-1/24*(kroneckerDelta(i,j)*kroneckerDelta(k,l)+kroneckerDelta(i,k)*kroneckerDelta(j,l)+kroneckerDelta(i,l)*kroneckerDelta(j,k))+1/6*(a2(i,j)*kroneckerDelta(k,l)+a2(i,k)*kroneckerDelta(j,l)+a2(i,l)*kroneckerDelta(j,k)+a2(k,l)*kroneckerDelta(i,j)+a2(j,l)*kroneckerDelta(i,k)+a2(j,k)*kroneckerDelta(i,l));
                end
            end
        end 
    end
end

