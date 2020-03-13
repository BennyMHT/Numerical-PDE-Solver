Node = zeros(n,n);
N = n*n;%number of nodes
Node(1:N) = 1:N;

f = zeros(625,1);
ng = (25-1)/6;%size of the grid for each block
for i = 1:4
    b1 = fix(source(i)./4)+1;
    b2 = rem(source(i),4);
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

force = reshape(f,25,25);

pcolor(force)