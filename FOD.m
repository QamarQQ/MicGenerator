function [ psi ,angle ] = untitled2( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i=1:2
    for j=1:2
        if i==j
            U(i,j)=1;
        else
            U(i,j)=0;
        end
    end
end

B=A-(1/2)*U;
p=cell(1,180);
f=cell(1,180);
f4=cell(1,180);

% A4=Cl_Natural(A);
% A4=Cl_Linear2D(A)
 A4=Cl_Quadratic(A);

%Second order
for theta=1:180
    p{theta}=[cosd(theta), sind(theta), 0];    
    for i=1:2
        for j=1:2
            %Second order f
            f{1,theta}(i,j)=p{theta}(i)*p{theta}(j)-(1/2)*U(i,j);
            for k=1:2
                for l=1:2
                    %4th order
                    %A4(i,j,k,l)=A(i,j)*A(k,l);
                    B4(i,j,k,l)=A4(i,j,k,l)-(1/6)*(U(i,j)*A(k,l)+U(i,k)*A(j,l)+...
                        U(i,l)*A(j,k)+U(j,k)*A(i,l)+U(j,l)*A(i,k)+U(k,l)*A(i,j))+...
                        (1/24)*(U(i,j)*U(k,l)+U(i,k)*U(j,l)+U(i,l)*U(j,k));
                    f4{1,theta}(i,j,k,l)=p{theta}(i)*p{theta}(j)*p{theta}(k)*p{theta}(l)-...
                        (1/6)*(U(i,j)*p{theta}(k)*p{theta}(l)+U(i,k)*p{theta}(j)*p{theta}(l)+...
                        U(i,l)*p{theta}(j)*p{theta}(k)+U(j,k)*p{theta}(i)*p{theta}(l)+...
                        U(j,l)*p{theta}(i)*p{theta}(k)+U(k,l)*p{theta}(i)*p{theta}(j))+...
                        (1/24)*(U(i,j)*U(k,l)+U(i,k)*U(j,l)+U(i,l)*U(j,k));
                    %psi
                    B4F4{theta}(i,j,k,l)=B4(i,j,k,l)*f4{theta}(i,j,k,l);
                    B4F4_1{theta}=sum(B4F4{theta});
                    B4F4_2{theta}=sum(B4F4_1{theta});
                    B4F4_3{theta}=sum(B4F4_2{theta});
                    B4F4_tot{theta}=sum(B4F4_3{theta});
%                     psi(theta)=1/(2*pi)+(2/pi)*(B(i,j)*f{theta}(i,j))+...
%                         (315/(32*pi))*(B(
                end
            end
        end
    end
    psi(theta)=1/(2*pi)+(2/pi)*...
        (B(1,1)*f{theta}(1,1)+B(1,2)*f{theta}(1,2)+...
        B(2,1)*f{theta}(2,1)+B(2,2)*f{theta}(2,2))+...
        (8/pi)*(B4F4_tot{theta});
    angle(theta)=deg2rad(theta);
end

% figure()
% hold on; box on; grid on
% plot(angle,psi)
% grid off
% xlabel ('$\theta$','Interpreter','latex')
% ylabel ('$\psi$','Interpreter','latex')
% axis([0,3.142,0,0.4])
% area=trapz(angle,psi)*2;
end

