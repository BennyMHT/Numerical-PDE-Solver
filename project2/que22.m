%specify parameters
roumax = 1;
vmax = 1;
d_x = 0.001;
d_t = 0.0005;
tao = 0.01;
theta0 = 0.5;
%initialize rou0 and v0 vector
N_x = 1/d_x+1;
rou0 = zeros(1,N_x);
v0 = zeros(1,N_x);
for i = 1:N_x
    x(i) = d_x*(i-1);
    if(x(i) <= 0.4) rou0(i) = 0.45;
    elseif(x(i) > 0.4 && x(i) <= 0.5) rou0(i) = 0.45+0.3*cos(5*pi*(x(i)-0.5));
    elseif(x(i) > 0.5 && x(i) <= 0.65) rou0(i) = 0.75;
    elseif(x(i) > 0.65 && x(i) <= 0.75) rou0(i) = 0.45+0.3*cos(5*pi*(x(i)-0.65));
    else rou0(i) = 0.45;
    end
end
for i = 1:N_x
    v0(i) = 1-rou0(i);
    %construct U(i+1/2),U(i-1/2),F(i+1/2),F(i-1/2),H vector
    U0(:,i) = [rou0(i);rou0(i)*v0(i)];
    F0(:,i) = [rou0(i)*v0(i);(v0(i)^2+theta0)*rou0(i)];
    H0(:,i) = [0;rou0(i)*(vmax*(1-rou0(i)/roumax)-v0(i))/tao];%ve = vmax(1-rou/roumax)
end
for i = 1:11
    [rou(i,:),v(i,:)] = laxsolver(U0,F0,H0,N_x,d_x,i*0.1,d_t,theta0,tao);
end

plot(x,rou(1,:),...
     x,rou(2,:),...
     x,rou(3,:),...
     x,rou(4,:),...
     x,rou(5,:),...
     x,rou(6,:),...
     x,rou(7,:),...
     x,rou(8,:),...
     x,rou(9,:),...
     x,rou(10,:),...
     x,rou(11,:))
legend('t = 0','t = 0.1','t = 0.2','t = 0.3','t = 0.4','t = 0.5','t = 0.6','t = 0.7','t = 0.8','t = 0.9','t = 1')
title('rou(x,t) theta0 = 0.5')

figure
plot(x,v(1,:),...
     x,v(2,:),...
     x,v(3,:),...
     x,v(4,:),...
     x,v(5,:),...
     x,v(6,:),...
     x,v(7,:),...
     x,v(8,:),...
     x,v(9,:),...
     x,v(10,:),...
     x,v(11,:))
legend('t = 0','t = 0.1','t = 0.2','t = 0.3','t = 0.4','t = 0.5','t = 0.6','t = 0.7','t = 0.8','t = 0.9','t = 1')
title('v(x,t) theta0 = 0.5')


%creat solver using Lax Scheme
function [rou,v] = laxsolver(U,F,H,N_x,d_x,t,d_t,theta0,tao)
N_t = t/d_t;
for n = 1:N_t
    %determine lamda,lamda = v+sqrt(theta0)
    for i = 1:N_x-1
        %v = U2/U1
        if(U(2,i)/U(1,i)+sqrt(theta0) >= U(2,i+1)/U(1,i+1)+sqrt(theta0))
            lamda(i) = U(2,i)/U(1,i)+sqrt(theta0);
        else lamda(i) = U(2,i+1)/U(1,i+1)+sqrt(theta0);
        end
    end
    %lamda(N_x) = lamda(1) by periodic
    %lax scheme
    for j = 2:N_x-1
        U(:,j) = U(:,j)-d_t/d_x*((F(:,j+1)-F(:,j-1))/2-lamda(j)*(U(:,j+1)-U(:,j))/2+lamda(j-1)*(U(:,j)-U(:,j-1))/2)...
                 +d_t*H(:,j);
        F(:,j) = [U(1,j)*U(2,j)/U(1,j);U(1,j)*((U(2,j)/U(1,j))^2+theta0)];%get F at n+1
        H(:,j) = [0;U(1,j)*(1-U(1,j)-U(2,j)/U(1,j))/tao];%get H at n+1
    end
end
rou = U(1,:);
v = U(2,:)./U(1,:);
end










