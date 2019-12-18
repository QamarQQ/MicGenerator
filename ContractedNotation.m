function [ a4_bar ] = ContractedNotation( A_bar )
%ContractedNotation Convert 4th order tensor to contracted notation
%   Convert 4th order tensor to contracted notation
    for m=1:6
        for n=1:6
            if (m==1 || m==2 || m==3) && (n==1 || n==2 || n==3)
                a4_bar(m,m,n,n)=A_bar(m,n);
            elseif (m==1 || m==2 || m==3) && n==4
                a4_bar(m,m,2,3)=A_bar(m,n);
                a4_bar(m,m,3,2)=A_bar(m,n);
            elseif (m==1 || m==2 || m==3) && n==5
                a4_bar(m,m,1,3)=A_bar(m,n);
                a4_bar(m,m,3,1)=A_bar(m,n);
            elseif (m==1 || m==2 || m==3) && n==6
                a4_bar(m,m,1,2)=A_bar(m,n);
                a4_bar(m,m,2,1)=A_bar(m,n);  
            elseif (n==1 || n==2 || n==3) && m==4
                a4_bar(2,3,n,n)=A_bar(m,n);
                a4_bar(3,2,n,n)=A_bar(m,n);
            elseif (n==1 || n==2 || n==3) && m==5
                a4_bar(1,3,n,n)=A_bar(m,n);
                a4_bar(3,1,n,n)=A_bar(m,n);
            elseif (n==1 || n==2 || n==3) && m==6
                a4_bar(1,2,n,n)=A_bar(m,n);
                a4_bar(2,1,n,n)=A_bar(m,n);
            elseif m==4 && n==4
                a4_bar(2,3,2,3)=A_bar(m,n);
                a4_bar(3,2,3,2)=A_bar(m,n);
                a4_bar(2,3,3,2)=A_bar(m,n);
                a4_bar(3,2,2,3)=A_bar(m,n);
            elseif m==5 && n==5
                a4_bar(1,3,1,3)=A_bar(m,n);
                a4_bar(3,1,3,1)=A_bar(m,n);
                a4_bar(1,3,3,1)=A_bar(m,n);
                a4_bar(3,1,1,3)=A_bar(m,n);
            elseif m==6 && n==6
                a4_bar(1,2,1,2)=A_bar(m,n);
                a4_bar(2,1,2,1)=A_bar(m,n);
                a4_bar(1,2,2,1)=A_bar(m,n);
                a4_bar(2,1,1,2)=A_bar(m,n);
            elseif m==4 && n==5
                a4_bar(2,3,1,3)=A_bar(m,n);
                a4_bar(3,2,3,1)=A_bar(m,n);
                a4_bar(2,3,3,1)=A_bar(m,n);
                a4_bar(3,2,1,3)=A_bar(m,n);
            elseif m==4 && n==6
                a4_bar(2,3,1,2)=A_bar(m,n);
                a4_bar(3,2,2,1)=A_bar(m,n);
                a4_bar(2,3,2,1)=A_bar(m,n);
                a4_bar(3,2,1,2)=A_bar(m,n);
            elseif m==5 && n==6
                a4_bar(1,3,1,2)=A_bar(m,n);
                a4_bar(3,1,2,1)=A_bar(m,n);
                a4_bar(1,3,2,1)=A_bar(m,n);
                a4_bar(3,1,1,2)=A_bar(m,n);
            elseif m==5 && n==4
                a4_bar(1,3,3,2)=A_bar(m,n);
                a4_bar(3,1,2,3)=A_bar(m,n);
                a4_bar(1,3,2,3)=A_bar(m,n);
                a4_bar(3,1,3,2)=A_bar(m,n);
            elseif m==6 && n==5
                a4_bar(2,1,1,3)=A_bar(m,n);
                a4_bar(1,2,3,1)=A_bar(m,n);
                a4_bar(2,1,3,1)=A_bar(m,n);
                a4_bar(1,2,1,3)=A_bar(m,n);
            elseif m==6 && n==4
                a4_bar(1,2,3,2)=A_bar(m,n);
                a4_bar(2,1,2,3)=A_bar(m,n);
                a4_bar(1,2,2,3)=A_bar(m,n);
                a4_bar(2,1,3,2)=A_bar(m,n);
            end
        end
    end
end
                

