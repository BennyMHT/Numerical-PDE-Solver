%specify parameters
B = [1,7,14,16];
n = 49;%size of the computational domain
hmin = 7;
w = 0.8;
v1 = 1;
v2 = 1;
vc = 1; 
t = 100;

%initialize vector phi
[A,f] = laplace(n,B);
phi = zeros(n*n,1);
for i = 1:t
    phi = vgridsolver(phi,f,n,B,hmin,v1,v2,vc,w);
    residue4(i) = norm(f-A*phi);
end

%convergence for 2 grid scheme 
for i = 1:t
    it_num(i) = (v1+v2)*(1+0.25^1+0.25^2+0.25^3)*i;
end

plot(it_num,log(residue4))
xlabel('times of iteration')
ylabel('log||ResV||')
title('convergence for 4 grid V cycle')

% solution = reshape(phi,n,n);
% figure
% [C,H] = contour(solution);
% clabel(C,H)
% title('Vgrid')
%%

%create vgrid solver
function  phi = vgridsolver(phi,f,n,B,hmin,v1,v2,vc,w)
%create matrix A on fine mesh h
A = laplace(n,B);

if~(n == hmin)
    %pre-smooth for v1 iterations
    phi = GS(phi,n*n,A,f,w,v1);
    r = f-A*phi;%compute residue on fine mesh
    [rc,nc] = restrict(reshape(r,n,n),n);
    ec = vgridsolver(zeros(nc*nc,1),reshape(rc,nc*nc,1),nc,B,hmin,v1,v2,vc,w);
    [e,n] = prolongate(reshape(ec,nc,nc),nc);
    phi = phi+reshape(e,n*n,1);
    %post-smooth for v2 iterations
    phi = GS(phi,n*n,A,f,w,v2);
else
    phi = GS(phi,n*n,A,f,w,vc);
end

end
%%

%function to create sparse laplacian operater matrix A and forcing
function [A,f] = laplace(n,B)
h = 1/(n-1);%length of the space unit
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
A = kron(D,I)+kron(I,D);
A = -A./h^2;

Node = zeros(n,n);
N = n*n;%number of nodes
Node(1:N) = 1:N;


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


%create forcing vector
f = zeros(N,1);
ng = (n-1)/6;%size of the grid for each block
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

end
%%

%creat Jacobi solver function
function phi1 = Jacobi(phi1,N,A,f,w,t)
Dia = diag(diag(A));
Up = -triu(A,1);
Lo = -tril(A,-1);
RJ = eye(N)-Dia^(-1)*A;%Jacobi iteration matrix
for i = 1:t
    phi1 = (w*RJ+(1-w)*eye(N))*phi1+w*Dia^(-1)*f;
end
end
%%

%create Gauss-Seidel solver function
function phi2 = GS(phi2,N,A,f,w,t)
Dia = diag(diag(A));
Up = -triu(A,1);
Lo = -tril(A,-1);
RGS = (Dia-Lo)^(-1)*Up;%Gauss-Seidel iteration matrix
for i =1:t
    phi2 = (w*RGS+(1-w)*eye(N))*phi2+w*(Dia-Lo)^(-1)*f;
end
end
%%

%restriction function
function [A,nc] = restrict(M,n)
nc = (n+1)/2;
for i = 1:nc
    temp1(i,:) = M(2*i-1,:);
end
for i = 1:nc
    A(:,i) = temp1(:,2*i-1);
end

end
%%

%prolongation function
function [A,n] = prolongate(M,nc)
n = 2*nc-1;
A = zeros(n,n);
%fill up 2i,2j
for i = 1:nc
    for j = 1:nc
    A(2*i-1,2*j-1) = M(i,j);
    end
end
%fill up 2i+1,2j
for i = 1:nc
    for j = 1:nc-1
        A(2*i-1,2*j) = (A(2*i-1,2*j-1)+A(2*i-1,2*j+1))/2;
    end
end
%fill up 2i,2j+1
for i = 1:nc-1
    for j = 1:nc
        A(2*i,2*j-1) = (A(2*i-1,2*j-1)+A(2*i+1,2*j-1))/2;
    end
end
%fill up 2i+1,2j+1
for i = 1:nc-1
    for j = 1:nc-1
        A(2*i,2*j) = (A(2*i-1,2*j-1)+A(2*i+1,2*j-1)+A(2*i-1,2*j+1)+A(2*i+1,2*j+1))/4;
    end
end

end









