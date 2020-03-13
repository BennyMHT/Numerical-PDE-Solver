function [solution_J,bot_j,top_j,lft_j,rit_j,solution_GS,bot_gs,top_gs,lft_gs,rit_gs] = iteration(n,B,w,t)
% n is size of the grid which has to be 6n+1,
% w is relaxaton factor which for this solver can only be chosen from [0,1]

%create sparse laplacian operater matrix A

ng = (n-1)/6;%size of the grid for each block
h = 1/(n-1);%length of the space unit
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
A = kron(D,I)+kron(I,D);
A = -A./h^2;

Node = zeros(n,n);
N = n*n;%number of nodes
Node(1:N) = 1:N;
%%

%specify boundary condition of matrix A

% ***bottom of the domain***
j = 1;
for i = 2:n-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i, j+1) ) = 0;
   A(ANode_i, Node(i+1, j) ) = 0;
   A(ANode_i, Node(i-1, j) ) = 0;
end

% ***top of the domain***
j = n;
for i = 2:n-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i, j-1) ) = 0;
   A(ANode_i, Node(i+1, j) ) = 0;
   A(ANode_i, Node(i-1, j) ) = 0;
end

% ***left of the domain***
i = 1;
for j = 2:n-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i, j+1) ) = 0;
   A(ANode_i, Node(i, j-1) ) = 0;
   A(ANode_i, Node(i+1, j) ) = 0;
end

% ***right of the domain***
i = n;
for j = 2:n-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i, j+1) ) = 0;
   A(ANode_i, Node(i, j-1) ) = 0;
   A(ANode_i, Node(i-1, j) ) = 0;
end

% domain corners condition

% ***bottom left point***
ANode_i = Node(1,1);
A(ANode_i, Node(2, 1) ) = 0;
A(ANode_i, Node(1, 2) ) = 0;

% ***bottom right point***
ANode_i = Node(1,n);
A(ANode_i, Node(2, n) ) = 0;
A(ANode_i, Node(1, n-1) ) = 0;

% ***top left point***
ANode_i = Node(n,1);
A(ANode_i, Node(n-1, 1) ) = 0;
A(ANode_i, Node(n, 2) ) = 0;

% ***top right point***
ANode_i = Node(n,n);
A(ANode_i, Node(n-1, n) ) = 0;
A(ANode_i, Node(n, n-1) ) = 0;


%derive D,U,L matrix from A Note:A = D-U-L
Dia = diag(diag(A));
Up = -triu(A,1);
Lo = -tril(A,-1);
%%

%create forcing matrix
f = zeros(N,1);
for i = 1:4
    b1 = fix(B(i)./4)+1;
    b2 = rem(B(i),4);
    if~(b2 == 0)
        j = 1+b1*ng;%row of the first point of the source block
        k = 1+b2*ng;%column of the first point
        f(Node(k:k+ng,j:j+ng),1) = 1;
    else
        j = 1+(b1-1)*ng;
        k = 1+4*ng;
        f(Node(k:k+ng,j:j+ng),1) = 1;
    end    
end
%%

%initialize vector phi
phi1 = zeros(N,1);
phi2 = zeros(N,1);

%creat Jacobi scheme
RJ = eye(N)-Dia^(-1)*A;%Jacobi iteration matrix
for i = 1:t
    e_j(i) = norm((w*RJ+(1-w)*eye(N))*phi1+w*Dia^(-1)*f-phi1);%get the error vector of Jocabi solver for each iteration
    phi1 = (w*RJ+(1-w)*eye(N))*phi1+w*Dia^(-1)*f;
end

%create Gauss-Seidel scheme
RGS = (Dia-Lo)^(-1)*Up;%Gauss-Seidel iteration matrix
for i =1:t
    e_gs(i) = norm((w*RGS+(1-w)*eye(N))*phi2+w*(Dia-Lo)^(-1)*f-phi2);%get the error vector of Jocabi solver for each iteration
    phi2 = (w*RGS+(1-w)*eye(N))*phi2+w*(Dia-Lo)^(-1)*f;
end

%reshape phi vector to matrix
solution_J = reshape(phi1,n,n);
solution_GS = reshape(phi2,n,n);

%create coordinate matrix
for i = 1:n
   for j = 1:n
       x(i,j)= (i-1)*h;
       y(i,j)= (j-1)*h;
   end
end

figure
[c,h] = contour(solution_J);
clabel(c,h)
%mesh(x,y,solution_J)
title('Jacobi Method')

figure
[c,h] = contour(solution_GS);
clabel(c,h)
%mesh(x,y,solution_GS)
title('Gauss Seidel Method')
%%

%show convergence for the solvers
it_num = 1:t;
%convergence for Jacobi solver
figure
plot(it_num,log(e_j))
xlabel('times of iteration')
ylabel('log||eJ||')
title('convergence fot Jacobi method')

%convergence for Gauss-Seidel solver
figure
plot(it_num,log(e_gs))
xlabel('times of iteration')
ylabel('log||eGS||')
title('convergence fot Gauss Seidel method')

%calculate flux on boundaries
%for Jacobi solver
bot_j = (-3*solution_J(1,:)+4*solution_J(2,:)-solution_J(3,:))/(2*h);
top_j = (-solution_J(n-2,:)+4*solution_J(n-1,:)-3*solution_J(n,:))/(2*h);
lft_j = (-3*solution_J(:,1)+4*solution_J(:,2)-solution_J(:,3))/(2*h);
rit_j = (-solution_J(:,n-2)+4*solution_J(:,n-1)-3*solution_J(:,n))/(2*h);

%for Gauss Seidel solver
bot_gs = (-3*solution_GS(1,:)+4*solution_GS(2,:)-solution_GS(3,:))/(2*h);
top_gs = (-solution_GS(n-2,:)+4*solution_GS(n-1,:)-3*solution_GS(n,:))/(2*h);
lft_gs = (-3*solution_GS(:,1)+4*solution_GS(:,2)-solution_GS(:,3))/(2*h);
rit_gs = (-solution_GS(:,n-2)+4*solution_GS(:,n-1)-3*solution_GS(:,n))/(2*h);


end












