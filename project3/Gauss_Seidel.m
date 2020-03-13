%==============question 2.b================
%convergence for Jacobi method
it_num = 1:5000;
for i = 1:9
    res = iteration(25,[1,7,14,16],0.2*i,5000);
    plot(it_num,log(res))
    hold on
end
xlabel('times of iteration')
ylabel('log||ResGS||')
title('convergence for Gauss Seidel method')
legend('w = 0.2','w = 0.4','w = 0.6','w = 0.8','w = 1','w = 1.2','w = 1.4','w = 1.6','w = 1.8')

it_num_GS = 1:5000;
res_GS = iteration(25,[1,7,14,16],0.8,5000);
%%

function [res_gs] = iteration(n,B,w,t)
% n is size of the grid which has to be 6n+1,
% w is relaxaton factor which for this solver can only be chosen from [0,2]

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
phi2 = zeros(N,1);

%create Gauss-Seidel scheme
RGS = (Dia-Lo)^(-1)*Up;%Gauss-Seidel iteration matrix
for i =1:t
    phi2 = (w*RGS+(1-w)*eye(N))*phi2+w*(Dia-Lo)^(-1)*f;
    res_gs(i) = norm(f - A*phi2);%get residue in each iteration
end

%reshape phi vector to matrix
solution = reshape(phi2,n,n);


end












