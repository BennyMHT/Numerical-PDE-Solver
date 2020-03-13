ve = 0.55;
t0 = 0.7;
c = 0.5;
t = 0:0.01:2*pi;
x = -(ve+sqrt(t0))+1+(ve+sqrt(t0))*cos(t);
y1 = -0.5*(ve+sqrt(t0))*sin(t);
y2 = -0.5*(ve-sqrt(t0))*sin(t);
plot(x,y1,x,y2)