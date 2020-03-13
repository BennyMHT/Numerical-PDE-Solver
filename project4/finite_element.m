% load grids 
load 'grids.mat'

% specify parameters
K = [0.4,0.6,0.8,1.1,1];
Bi = 0.1;

% read numbers of elements from different mesh
for i = 1:7
    [M,N] = size(coarse.theta{i});
    coa_nofe(i) = M;
    [M,N] = size(medium.theta{i});
    med_nofe(i) = M;
    [M,N] = size(fine.theta{i});
    fin_nofe(i) = M;
end

% ***** medium mesh case *****

% initialization
Num = medium.nodes;
A = sparse(Num,Num);
F = sparse(Num,1);
Alocal = zeros(3);
Alocal_ex = zeros(2);
Flocal_root = zeros(2,1);

% fill out matrix of interior domain
for i = 1:5
    for k = 1:med_nofe(i)
        n1 = medium.theta{i}(k,1);
        n2 = medium.theta{i}(k,2);
        n3 = medium.theta{i}(k,3);
        x1 = medium.coor(n1,1);
        y1 = medium.coor(n1,2);
        x2 = medium.coor(n2,1);
        y2 = medium.coor(n2,2);
        x3 = medium.coor(n3,1);
        y3 = medium.coor(n3,2);
        coef = [1 x1 y1;1 x2 y2;1 x3 y3];
        rhs = [1;0;0];
        c1 = coef\rhs;
        cx1 = c1(2);
        cy1 = c1(3);
        rhs = [0;1;0];
        c2 = coef\rhs;
        cx2 = c2(2);
        cy2 = c2(3);
        rhs = [0;0;1];
        c3 = coef\rhs;
        cx3 = c3(2);
        cy3 = c3(3);
        Alocal(1,1) = cx1*cx1+cy1*cy1;
        Alocal(1,2) = cx1*cx2+cy1*cy2;
        Alocal(1,3) = cx1*cx3+cy1*cy3;
        Alocal(2,1) = cx2*cx1+cy2*cy1;
        Alocal(2,2) = cx2*cx2+cy2*cy2;
        Alocal(2,3) = cx2*cx3+cy2*cy3;
        Alocal(3,1) = cx3*cx1+cy3*cy1;
        Alocal(3,2) = cx3*cx2+cy3*cy2;
        Alocal(3,3) = cx3*cx3+cy3*cy3;
        area = abs(det(coef)/2);
        Alocal = K(i)*area*Alocal;
        for a = 1:3
            for b = 1:3
                I = medium.theta{i}(k,a);
                J = medium.theta{i}(k,b);
                A(I,J) = A(I,J)+Alocal(a,b);
            end
        end
    end
end

% fill out boundary with robin bc
for k = 1:med_nofe(6)
    n1 = medium.theta{6}(k,1);
    n2 = medium.theta{6}(k,2);
    x1 = medium.coor(n1,1);
    y1 = medium.coor(n1,2);
    x2 = medium.coor(n2,1);
    y2 = medium.coor(n2,2);
    h = sqrt((x1-x2)^2+(y1-y2)^2);
    Alocal_ex = [1/3 1/6;1/6 1/3];
    Alocal_ex = Bi*h*Alocal_ex;
    for a = 1:2
        for b = 1:2
            I = medium.theta{6}(k,a);
            J = medium.theta{6}(k,b);
            A(I,J) = A(I,J)+Alocal_ex(a,b);
        end
    end
end

% fill out load vector
for k = 1:med_nofe(7)
    n1 = medium.theta{7}(k,1);
    n2 = medium.theta{7}(k,2);
    x1 = medium.coor(n1,1);
    y1 = medium.coor(n1,2);
    x2 = medium.coor(n2,1);
    y2 = medium.coor(n2,2);
    h = sqrt((x1-x2)^2+(y1-y2)^2);
    Flocal_root = [1/2;1/2];
    Flocal_root = h*Flocal_root;
    for a = 1:2
        I = medium.theta{7}(k,a);
        F(I) = F(I)+Flocal_root(a);
    end
end

% solve for u
u = A\F;
% get mean temperature at root
T = F'*u;

plotsolution(medium,u)





