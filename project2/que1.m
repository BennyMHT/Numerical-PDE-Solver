%solve and plot velocity vector in each initial condition
%creat velocity vector
N_j = 200;
u0 = zeros(1,N_j);

%case a single soliton u(x,0) = u1(x,0) 
%initiate velocity vector 
for j = 1:N_j
    x(j) = -10+(j-1)*0.1;
    u0(j) = -20/(2*(cosh(sqrt(20)*x(j)/2))^2);
end
%plot initial velocity vector
figure
subplot(2,1,1)
plot(x,u0)
title('case a')

%solve velocity vector by step
ua = solitonsolver(N_j,2,0.001,u0);
%plot final velocity vector
subplot(2,1,2)
plot(x,ua)

%case b single soliton u(x,0) = -10*exp(-x^2) 
%initiate velocity vector 
for j = 1:N_j
    u0(j) = -10*exp(-x(j)^2);
end
%plot initial velocity vector
figure
subplot(2,1,1)
plot(x,u0)
title('case b')

%solve velocity vector by step
ub = solitonsolver(N_j,2,0.001,u0);
%plot final velocity vector
subplot(2,1,2)
plot(x,ub)

%case c two_soliton soliton u(x,0) = -6/(cosh(x)).^2 
%initiate velocity vector 
for j = 1:N_j
    u0(j) = -6/(cosh(x(j)))^2;
end
%plot initial velocity vector
figure
subplot(2,1,1)
plot(x,u0)
title('case c')

%solve velocity vector by step
uc = solitonsolver(N_j,2,0.001,u0);
%plot final velocity vector
subplot(2,1,2)
plot(x,uc)

%case d own two_soliton soliton where v1=14,v2=6,both x0=0
%u(x,0) = -14/(2*(cosh(sqrt(14)*x/2))^2)-6/(2*(cosh(sqrt(6)*x/2))^2)
%initiate velocity vector 
for j = 1:N_j
    u0(j) = -14/(2*(cosh(sqrt(14)*x(j)/2))^2)-6/(2*(cosh(sqrt(6)*x(j)/2))^2);
end
%plot initial velocity vector
figure
subplot(2,1,1)
plot(x,u0)
title('case d')

%solve velocity vector by step
ud = solitonsolver(N_j,2,0.001,u0);
%plot final velocity vector
subplot(2,1,2)
plot(x,ud)

%case e own two_soliton soliton where (v1,x0)=(14,-3),(v2,x0)=(6,3)
%u(x,0) = -14/(2*(cosh(sqrt(14)*(x+3)/2))^2)-6/(2*(cosh(sqrt(6)*(x-3)/2))^2)
%initiate velocity vector 
for j = 1:N_j
    u0(j) = -14/(2*(cosh(sqrt(14)*(x(j)+3)/2))^2)-6/(2*(cosh(sqrt(6)*(x(j)-3)/2))^2);
end
%plot initial velocity vector
figure
subplot(5,1,1)
plot(x,u0)
title('case e')

%solve velocity vector by step 
for i = 1:4
    ue(i,:) = solitonsolver(N_j,i*0.5,0.001,u0);
end
%plot final velocity vector
subplot(5,1,2)
plot(x,ue(1,:))
subplot(5,1,3)
plot(x,ue(2,:))
subplot(5,1,4)
plot(x,ue(3,:))
subplot(5,1,5)
plot(x,ue(4,:))

%case f u0 = 2*sin(x/pi)
%initiate velocity vector 
for j = 1:N_j
    u0(j) = 2*sin(x(j)/pi);
end
%plot initial velocity vector
figure
subplot(2,1,1)
plot(x,u0)
title('case f')

%solve velocity vector by step
uf = solitonsolver(N_j,2,0.001,u0);
%plot final velocity vector
subplot(2,1,2)
plot(x,uf)


%solver of the governing equation
function [solution] = solitonsolver(J,t,dt,u)
%calculate number of time steps
N_t = t/dt;

%creat spatial discretize function f
function [sp] = f(un)
for j = 1:J
    if(j == 1)
        sp(j) = -(-un(J-1)+2*un(J)-2*un(j+1)+un(j+2))/(2*0.1^3)+6*un(j)*(un(j+1)-un(J))/(2*0.1);
    elseif(j == 2)
        sp(j) = -(-un(J)+2*un(1)-2*un(j+1)+un(j+2))/(2*0.1^3)+6*un(j)*(un(j+1)-un(j-1))/(2*0.1);
    elseif(j == J-1)
        sp(j) = -(-un(j-2)+2*un(j-1)-2*un(j+1)+un(1))/(2*0.1^3)+6*un(j)*(un(j+1)-un(j-1))/(2*0.1);
    elseif(j == J)
        sp(j) = -(-un(j-2)+2*un(j-1)-2*un(1)+un(2))/(2*0.1^3)+6*un(j)*(un(1)-un(j-1))/(2*0.1);
    else
        sp(j) = -(-un(j-2)+2*un(j-1)-2*un(j+1)+un(j+2))/(2*0.1^3)+6*un(j)*(un(j+1)-un(j-1))/(2*0.1);
    end
end
end

%time integration using Runge-Kutta scheme
for n = 1:N_t
    a1 = dt*f(u);
    a2 = dt*f(u+a1/2);
    a3 = dt*f(u+a2/2);
    a4 = dt*f(u+a3);
    u = u+(a1+a2+a3+a4)/6;    
end

solution = u;%return velocity field

end %end solver











