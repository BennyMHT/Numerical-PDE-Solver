%specify parameters
b = 1:16;
B = nchoosek(b,4);
n = 25;%size of the computational domain
w = 0.8;
v1 = 1;
v2 = 1;
vc = 4; 
t = 500;

%plot convergence of multigrid scheme
res_M4 = multigridsolver([1,7,14,16],n,w,v1,v2,vc,t);
for i = 1:t
    it_num_M(i) = (v1+v2+vc/4)*i;
end
% plot(it_num_M,log(res_M(1,:)),it_num_M,log(res_M(2,:)))
% xlabel('times of iteration')
% ylabel('log||ResM||')
% title('convergence for Multigrid Gauss Seidel method')
% legend('w = 0.5','w = 0.8')
%%

%creat loop to calculate force location given unknow boundary flux
% flux = load('flux_unknown.mat');
% flux=flux.flux;
% num = nchoosek(16,4);
% time = 0;%number of loop
% for i = 1:num
%     F = multigridsolver(B(i,:),n,w,v1,v2,t);
%     if (abs(F(10,1)-flux(10,1)) < 1e-4 && ...
%         abs(F(11,2)-flux(11,2)) < 1e-4 && ...
%         abs(F(12,3)-flux(12,3)) < 1e-4 && ...
%         abs(F(13,4)-flux(13,4)) < 1e-4)
%         source = B(i,:);
%         time = time+1;
%         break;
%     else
%         time = time+1;
%     end
% end

%%

%create multigrid solver
function residue = multigridsolver(B,n,w,v1,v2,vc,t)
%create matrix A on fine mesh h
[A,f] = laplace(n,B);

%create A matrix on a coarse grid 2h
nc = (n+1)/2;%size of the coarse grid
[A2,f2] = laplace(nc,B);

%initialize vector phi
phi = zeros(n*n,1);

%implement multigrid method
for time = 1:t
%pre-smooth for v1 iterations
[phi,res] = GS(phi,n*n,A,f,w,v1);
r = f-A*phi;%compute residue on fine mesh

%restriction
R = reshape(r,n,n);
for i = 1:nc
    temp1(i,:) = R(2*i-1,:);
end
for i = 1:nc
    R2(:,i) = temp1(:,2*i-1);
end
r2 = reshape(R2,nc*nc,1);

%solve for error on coarse grid
e2 = zeros(nc*nc,1);
[e2,res_e] = GS(e2,nc*nc,A2,r2,w,vc);

%prolongation
E2 = reshape(e2,nc,nc);
E = zeros(n,n);

% for j = 1:nc-1
%     for i = 1:nc-1
%         E(2*i,2*j) = E2(i,j);
%         E(2*i+1,2*j) = 0.5*(E2(i,j)+E2(i+1,j));
%         E(2*i,2*j+1) = 0.5*(E2(i,j)+E2(i,j+1));
%         E(2*i+1,2*j+1) = 0.25*(E2(i,j)+E2(i+1,j)+E2(i,j+1)+E2(i+1,j+1));
%     end
% end


%fill up 2i,2j
for i = 1:nc
    for j = 1:nc
    E(2*i-1,2*j-1) = E2(i,j);
    end
end
%fill up 2i+1,2j
for i = 1:nc
    for j = 1:nc-1
        E(2*i-1,2*j) = (E(2*i-1,2*j-1)+E(2*i-1,2*j+1))/2;
    end
end
%fill up 2i,2j+1
for i = 1:nc-1
    for j = 1:nc
        E(2*i,2*j-1) = (E(2*i-1,2*j-1)+E(2*i+1,2*j-1))/2;
    end
end
%fill up 2i+1,2j+1
for i = 1:nc-1
    for j = 1:nc-1
        E(2*i,2*j) = (E(2*i-1,2*j-1)+E(2*i+1,2*j-1)+E(2*i-1,2*j+1)+E(2*i+1,2*j+1))/4;
    end
end

e = reshape(E,n*n,1);

phi = phi+e;

%post-smooth for v2 iterations
[phi,residue(time)] = GS(phi,n*n,A,f,w,v2);%get residue for each cycle
solution = reshape(phi,n,n);
end

% %plot contour
% figure
% [C,H] = contour(solution);
% clabel(C,H)
% title('Multigrid Gauss Seidel Method')

% %calculate flux on boundaries
% h = 1/(n-1);
% bot = (-3*solution(1,:)+4*solution(2,:)-solution(3,:))/(2*h);
% top = (-solution(n-2,:)+4*solution(n-1,:)-3*solution(n,:))/(2*h);
% lft = (-3*solution(:,1)+4*solution(:,2)-solution(:,3))/(2*h);
% rit = (-solution(:,n-2)+4*solution(:,n-1)-3*solution(:,n))/(2*h);
% %form flux matrix
% flux = zeros(n,4);
% flux(:,1) = lft;
% flux(:,2) = rit;
% flux(:,3) = bot';
% flux(:,4) = top';

end
%%

%function to create sparse laplacian operater matrix A and foring vector
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
function [phi1,res_j] = Jacobi(phi1,N,A,f,w,t)
Dia = diag(diag(A));
Up = -triu(A,1);
Lo = -tril(A,-1);
RJ = eye(N)-Dia^(-1)*A;%Jacobi iteration matrix
for i = 1:t
    phi1 = (w*RJ+(1-w)*eye(N))*phi1+w*Dia^(-1)*f;
    res_j = norm(f - A*phi1);%get residue in each iteration
end
end

%create Gauss-Seidel solver function
function [phi2,res_gs] = GS(phi2,N,A,f,w,t)
Dia = diag(diag(A));
Up = -triu(A,1);
Lo = -tril(A,-1);
RGS = (Dia-Lo)^(-1)*Up;%Gauss-Seidel iteration matrix
for i =1:t
    phi2 = (w*RGS+(1-w)*eye(N))*phi2+w*(Dia-Lo)^(-1)*f;
    res_gs = norm(f-A*phi2);%get residue for each iteration
end
end









