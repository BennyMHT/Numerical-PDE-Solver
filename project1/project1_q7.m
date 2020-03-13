i = 0;
for h = 0.1:0.1:2
    for b = 0.1:0.1:2
        if h+b <= 3.25
        i = i+1;
        [Q(i),I(i)] = Q_I_output(21,21,6.5,b,h);
        end
    end
end
k = convhull(Q,I);
plot(Q(k),I(k),Q,I,'.');