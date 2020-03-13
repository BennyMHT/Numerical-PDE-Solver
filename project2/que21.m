%parameters specification
roumax = 1;
vmax = 1;
d_x = 0.001;
d_t = 0.8*d_x/vmax;
%initiate density vector
%space domain [0,1]
N_x = 1/0.001+1;
rou0 = zeros(1,N_x);
for i = 1:N_x
    x(i) = d_x*(i-1);
    if(x(i) <= 0.4) rou0(i) = 0.45;
    elseif(x(i) > 0.4 && x(i) <= 0.5) rou0(i) = 0.45+0.3*cos(5*pi*(x(i)-0.5));
    elseif(x(i) > 0.5 && x(i) <= 0.65) rou0(i) = 0.75;
    elseif(x(i) > 0.65 && x(i) <= 0.75) rou0(i) = 0.45+0.3*cos(5*pi*(x(i)-0.65));
    else rou0(i) = 0.45;
    end
end

%solve continuity equation
for j = 1:11
rout(j,:) = trafficflow(rou0,N_x,d_x,0.1*(j-1),d_t,roumax,vmax);
end
plot(x,rout(1,:),...
     x,rout(2,:),...
     x,rout(3,:),...
     x,rout(4,:),...
     x,rout(5,:),...
     x,rout(6,:),...
     x,rout(7,:),...
     x,rout(8,:),...
     x,rout(9,:),...
     x,rout(10,:),...
     x,rout(11,:))
legend('t = 0','t = 0.1','t = 0.2','t = 0.3','t = 0.4','t = 0.5','t = 0.6','t = 0.7','t = 0.8','t = 0.9','t = 1')

%creat solver for finite volume scheme of traffic flow
function [rou] = trafficflow(rou,N_x,d_x,t,d_t,roumax,vmax)
N_t = t/d_t;
for n = 1:N_t
    for i = 2:N_x-1
        if(rou(i+1)+rou(i) >= roumax)
            rou(i) = rou(i)-d_t*(rou(i+1)*vmax*(1-rou(i+1)/roumax)-rou(i)*vmax*(1-rou(i)/roumax))/d_x;%f(rou) = rou*vmax*(1-rou/roumax)
        elseif(rou(i+1)+rou(i) <= roumax)
            rou(i) = rou(i)-d_t*(rou(i)*vmax*(1-rou(i)/roumax)-rou(i-1)*vmax*(1-rou(i-1)/roumax))/d_x;
        end
    end
end
end



