% case bb = 0.5
[q81_1,sol81_1] = Q_Sol_output(81,81,6.5,0.5,2);
[q41_1,sol41_1] = Q_Sol_output(41,41,6.5,0.5,2);
[q21_1,sol21_1] = Q_Sol_output(21,21,6.5,0.5,2);
[q11_1,sol11_1] = Q_Sol_output(11,11,6.5,0.5,2);

% restrict the solution obtained at the reference mesh to the current coarse mesh
for i = 1:11
    for j = 1:11
        e11_1(i,j) = abs(sol11_1(i,j)-sol81_1(1+(i-1)*8,1+(j-1)*8));
    end
end

for i = 1:21
    for j = 1:21
        e21_1(i,j) = abs(sol21_1(i,j)-sol81_1(1+(i-1)*4,1+(j-1)*4));
    end
end

for i = 1:41
    for j = 1:41
        e41_1(i,j) = abs(sol41_1(i,j)-sol81_1(1+(i-1)*2,1+(j-1)*2));
    end
end

% reshape solution matrix into vector
e11_1 = e11_1(:);
e12_1 = e21_1(:);
e41_1 = e41_1(:);

% calculate the errors of the flowrate
eq_11_1 = abs(q11_1-q81_1);
eq_21_1 = abs(q21_1-q81_1);
eq_41_1 = abs(q41_1-q81_1);

eq1 = [eq_11_1,eq_21_1,eq_41_1];
xeq1 = [1/10,1/20,1/40];
a1_1 = polyfit(log10(xeq1),log10(eq1),1);

% calculate the errors of the L2 norm

eL2_11_1 = 1/10*norm(e11_1,2);
eL2_21_1 = 1/20*norm(e21_1,2);
eL2_41_1 = 1/40*norm(e41_1,2);

eL2_1 = [eL2_11_1,eL2_21_1,eL2_41_1];
xeL2_1 = [1/10,1/20,1/40];
a2_1 = polyfit(log10(xeL2_1),log10(eL2_1),1);

% calculate the errors of the Linf norm
eLinf_11_1 = norm(e11_1,inf);
eLinf_21_1 = norm(e21_1,inf);
eLinf_41_1 = norm(e41_1,inf);

eLinf_1 = [eLinf_11_1,eLinf_21_1,eLinf_41_1];
xeLinf_1 = [1/10,1/20,1/40];
a3_1 = polyfit(log10(xeLinf_1),log10(eLinf_1),1);

% case bb = 0.7
[q81_2,sol81_2] = Q_Sol_output(81,81,6.5,0.7,2);
[q41_2,sol41_2] = Q_Sol_output(41,41,6.5,0.7,2);
[q21_2,sol21_2] = Q_Sol_output(21,21,6.5,0.7,2);
[q11_2,sol11_2] = Q_Sol_output(11,11,6.5,0.7,2);

% restrict the solution obtained at the reference mesh to the current coarse mesh
for i = 1:11
    for j = 1:11
        e11_2(i,j) = abs(sol11_2(i,j)-sol81_2(1+(i-1)*8,1+(j-1)*8));
    end
end

for i = 1:21
    for j = 1:21
        e21_2(i,j) = abs(sol21_2(i,j)-sol81_2(1+(i-1)*4,1+(j-1)*4));
    end
end

for i = 1:41
    for j = 1:41
        e41_2(i,j) = abs(sol41_2(i,j)-sol81_2(1+(i-1)*2,1+(j-1)*2));
    end
end

% reshape solution matrix into vector
e11_2 = e11_2(:);
e12_2 = e21_2(:);
e41_2 = e41_2(:);

% calculate the errors of the flowrate
eq_11_2 = abs(q11_2-q81_2);
eq_21_2 = abs(q21_2-q81_2);
eq_41_2 = abs(q41_2-q81_2);

eq2 = [eq_11_2,eq_21_2,eq_41_2];
xeq2 = [1/10,1/20,1/40];
a1_2 = polyfit(log10(xeq2),log10(eq2),1);

% calculate the errors of the L2 norm

eL2_11_2 = 1/10*norm(e11_2,2);
eL2_21_2 = 1/20*norm(e21_2,2);
eL2_41_2 = 1/40*norm(e41_2,2);

eL2_2 = [eL2_11_2,eL2_21_2,eL2_41_2];
xeL2_2 = [1/10,1/20,1/40];
a2_2 = polyfit(log10(xeL2_2),log10(eL2_2),1);

% calculate the errors of the Linf norm
eLinf_11_2 = norm(e11_2,inf);
eLinf_21_2 = norm(e21_2,inf);
eLinf_41_2 = norm(e41_2,inf);

eLinf_2 = [eLinf_11_2,eLinf_21_2,eLinf_41_2];
xeLinf_2 = [1/10,1/20,1/40];
a3_2 = polyfit(log10(xeLinf_2),log10(eLinf_2),1);

% plot
subplot(2,3,1)
loglog(xeq1,eq1);
title('q, b = 0.5')
grid on

subplot(2,3,2)
loglog(xeL2_1,eL2_1);
title('L2, b = 0.5')
grid on

subplot(2,3,3)
loglog(xeLinf_1,eLinf_1);
title('Linf, b = 0.5')
grid on

subplot(2,3,4)
loglog(xeq2,eq2);
title('q, b = 0.7')
grid on

subplot(2,3,5)
loglog(xeL2_2,eL2_2);
title('L2, b = 0.7')
grid on

subplot(2,3,6)
loglog(xeLinf_2,eLinf_2);
title('Linf, b = 0.7')
grid on




