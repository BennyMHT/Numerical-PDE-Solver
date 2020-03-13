function [FlowRate, I_xx] = ChannelFlow(N_xi, N_eta, ll, bb, hh)
%==================================================================================
% ChannelFlow.m
%==================================================================================
% The following code is the solution skeleton to the channel flow problem.
% Since this code is intended as a skeleton and should be easy to read,
% we chose not to implement the code using vectorization.
% This will cause slightly lower performance.
%
% Inputs:
% =======
%    N_xi  : Number Of nodes in the xi direction.
%    N_eta : Number of nodes int he eta direction
%    ll    : Channel perimeter
%    bb    : Width of the base of the channel
%    hh    : Height of the channel
%
% Outputs:
% ========
%    Solution : A xi-eta matrix containing the solution
%    Flowrate : The flowrate throught the channel
%    I_xx     : The moment of inertia of the channel
%
% Some additional information for the function:
%
% In this pseudo code, the following conventions are used:
%
%    1) The computational domain xi-eta has values of (0,0) corresponding to the
%       bottom left corner of the physical and computational domains. Correspondingly,
%       the upper-right corner in the xi-eta domain has a (xi, eta) value of (1,1).
%
%    2) We use a node map in the matrix 'Node'. Node(i,j) grabs node reference
%       number. Using the Node matrix, we can easily stamp/stencil the related
%       values into the finite difference matrix. To determine how the Node matrix
%       looks, you can easily type the creation commands in the command prompt.
%
% Originally written by: D.J.Willis
% Modified and distributed by A. Uranga with permission
%
%
%==================================================================================
%close all;


%%=== Geometric Parameters ========================================================

cc = bb+((ll-2*bb)^2/4-hh^2)^0.5;



%%=== Grid Details and parameters =================================================

d_xi  = 1./(N_xi-1);
d_eta = 1./(N_eta-1);
NumNodes = N_xi * N_eta;


%%=== Initializations =============================================================

%- Initializing the sparse matrix 'A'
%  Note: It is essential to allocate a sparse A-matrix due to memory restrictions.
A   = spalloc(NumNodes, NumNodes, 9*NumNodes);

%- Initializing the RHS.
RHS = ones(NumNodes,1);


%%=== Note numbering scheme =======================================================
%
% Creation of a node numbering scheme. The node numbering scheme is created here
% in order to simplify the overall solution process. The idea is as follows:
%
% We construct a matrix called Node, which has elements corresponding to the node
% numbers in the grid representation of the solution domain. This allows us to cycle
% through the resulting matrix, grab element (i,j) and easily find the (i +/- 1),
% and (j +/- 1) node numbers.
%
Node = zeros(N_xi, N_eta);
Node(1:NumNodes) = 1:NumNodes;


%%=== Jacobian ====================================================================
'Constructing The Jacobian'

for i = 1:N_xi
   for j = 1:N_eta
      xi(i,j)  = (i-1)*d_xi;
      eta(i,j) = (j-1)*d_eta;
      J(i,j)   = hh*((cc-bb)*eta(i,j)+bb);
   end
end

%%=== "A" Matrix ==================================================================
'Constructing the "A" Matrix'

%----------------------------------------------------------------------------------
%------------------------ INNER REGION OF THE DOMAIN ------------------------------
% We begin by considering the Inner region of the domain.
% This is essentially the fill-in for all nodes not touching the boundary.
% The boundary nodes are handled later.

for i = 2:N_xi-1
   for j = 2:N_eta-1
      ANode_i = Node(i,j);            % Setting A-Matrix position for node i,j
		%----------------------------------------------------------------------------
		%--------------------- The Transformation -----------------------------------
      % The various components of the transformation

		a = (cc-bb)^2*xi(i,j)^2+hh^2;
		b = ((cc-bb)*eta(i,j)+bb)*(cc-bb)*xi(i,j);
		c = ((cc-bb)*eta(i,j)+bb)^2;
		alpha = -2*b*(cc-bb);
		beta  = 0;
		d = 0;
		e = -alpha*hh/J(i,j);

		%----------------------------------------------------------------------------
		%--------------------- FILLING UP THE A MATRIX ------------------------------
		% The filling of the matrix is done via the stamping of the computational
      % molecule in the appropriate parts of the A-Matrix.


		%----------------------------------------------------------------------------
		%---------------------  RHS Part Of Computational Molecule  -----------------
		% Here A(p,q) is such that:
      %    p = the current node number on the grid, at which the
		%        differential equation is being approximated
      %    q = the neighboring point to the current node.
      % So, here we are using a stamping stencil based on the ANode_i matrix.
      % It may be worthwhile to display a reduced dimension version of ANode_i
      % to fully grasp what is happening here.

      A(ANode_i, Node(i+1, j+1) ) = -(-2*b/(4*d_xi*d_eta))/J(i,j)^2;
      A(ANode_i, Node(i+1, j  ) ) = -(e/(2*d_xi)+a/d_xi^2)/J(i,j)^2;
	  A(ANode_i, Node(i+1, j-1) ) = -(2*b/(4*d_xi*d_eta))/J(i,j)^2;


      %----------------------------------------------------------------------------
      %--------------------  Middle Part Of Computational Molecule  ---------------

      A(ANode_i, Node(i  , j+1) ) = -(d/(2*d_eta)+c/d_eta^2)/J(i,j)^2;
      A(ANode_i, Node(i  , j  ) ) = -(-2*a/d_xi^2-2*c/d_eta^2)/J(i,j)^2;
      A(ANode_i, Node(i  , j-1) ) = -(-d/(2*d_eta)+c/d_eta^2)/J(i,j)^2;


      %----------------------------------------------------------------------------
      %---------------------  LHS  Part Of Computational Molecule  ----------------

      A(ANode_i, Node(i-1, j+1) ) = -(2*b/(4*d_xi*d_eta))/J(i,j)^2;
	  A(ANode_i, Node(i-1, j  ) ) = -(-e/(2*d_xi)+a/d_xi^2)/J(i,j)^2;
	  A(ANode_i, Node(i-1, j-1) ) = -(-2*b/(4*d_xi*d_eta))/J(i,j)^2;

	end
end

%----------------------------------------------------------------------------------
%------------------------ BOUNDARY CONDITIONS -------------------------------------

%------------------------ Bottom of the domain ------------------------------------
j = 1;
for i = 2:N_xi-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i,j)) = 1;
   RHS(ANode_i) = 0;
end

%------------------------ Top of the domain ---------------------------------------
j = N_eta;
for i = 2:N_xi-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i, j-2) ) = 1/(2*d_eta)*((cc-bb)*eta(i,j)+bb)/J(i,j);
   A(ANode_i, Node(i, j-1) ) = -4/(2*d_eta)*((cc-bb)*eta(i,j)+bb)/J(i,j);
   A(ANode_i, Node(i, j  ) ) = 3/(2*d_eta)*((cc-bb)*eta(i,j)+bb)/J(i,j);
   A(ANode_i, Node(i-1, j) ) = 1/(2*d_xi)*(cc-bb)*xi(i,j)/J(i,j);
   A(ANode_i, Node(i+1, j) ) = -1/(2*d_xi)*(cc-bb)*xi(i,j)/J(i,j);
   RHS(ANode_i) = 0;
   % Fill in the boundary condition for the top of the domain

end

%------------------------ Left side of the domain ---------------------------------
i = 1;
for j = 2:N_eta-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i  , j) ) = -3/(2*d_xi)*hh/J(i,j);
   A(ANode_i, Node(i+1, j) ) = 4/(2*d_xi)*hh/J(i,j);
   A(ANode_i, Node(i+2, j) ) = -1/(2*d_xi)*hh/J(i,j);
   RHS(ANode_i) = 0;
   % Fill in the boundary condition for the left of the domain

end

%------------------------ Right side of the domain --------------------------------
i = N_xi;
for j = 2:N_eta-1
   ANode_i = Node(i,j);
   A(ANode_i, Node(i,j)) = 1;
   RHS(ANode_i) = 0;
   % Fill in the boundary condition for the Right side of the domain

end

%------------------------ DOMAIN CORNERS ------------------------------------------

%------------------------ Bottom left ---------------------------------------------
ANode_i = Node(1,1);
A(ANode_i, Node(1,1)) = 1;
RHS(ANode_i) = 0;
% Fill in the boundary condition for the bottom left corner of the domain


%------------------------ Bottom right --------------------------------------------
ANode_i = Node(N_xi,1);
A(ANode_i, Node(N_xi,1)) = 1;
RHS(ANode_i) = 0;
% Fill in the boundary condition for the bottom right corner of the domain


%------------------------ Top Left ------------------------------------------------
ANode_i = Node(1,N_eta);  % Setting A_Matrix position for node i,j
A(ANode_i, Node(1  , N_eta) ) = -3/(2*d_xi)*hh/J(1,N_eta);
A(ANode_i, Node(1+1, N_eta) ) = 4/(2*d_xi)*hh/J(1,N_eta);
A(ANode_i, Node(1+2, N_eta) ) = -1/(2*d_xi)*hh/J(1,N_eta);
RHS(ANode_i) = 0;
% Fill in the boundary condition for the top left corner of the domain


%------------------------ Top right -----------------------------------------------
ANode_i = Node(N_xi,N_eta);
A(ANode_i, Node(N_xi,N_eta)) = 1;
RHS(ANode_i) = 0;
% Fill in the boundary condition for the top right corner of the domain


%%=== Solving the system Ax=b =====================================================
'Solving the system'

Sol = A\RHS;
Solution = reshape(Sol,N_xi,N_eta);


%%=== Post-processing the solution ================================================
'Post-processing'

%----------------------------------------------------------------------------------
%------------------------ Computing the flow rate & moment of inertia -------------
% With the velocity known, we can compute the flowrate and the moment of inertia.
% As a check of your flow rate integral try computing the area of the channel
% and compare that with the analytical result.

%-- Flow rate
FlowRate = 0;
for i = 1:N_xi-1
    for j = 1:N_eta-1
        FlowRate = FlowRate+(Solution(i,j)+Solution(i+1,j)+Solution(i,j+1)+Solution(i+1,j+1))/4*...
        (J(i,j)+J(i+1,j)+J(i,j+1)+J(i+1,j+1))/4*d_xi*d_eta;
    end
end
FlowRate = FlowRate*2;
% Fill In Your Computations For the Flow Rate


%-- Moment of inertia
t = 0.05;
I_xx = bb*t^3/6+bb*t*hh^2/2+t*((cc-bb)^2+hh^2)^0.5/6*(hh^2+t^2*(cc-bb)^2/((cc-bb)^2+hh^2));
% Fill In Your Computations For the Moment Of Inertia



%%=== Information for the plots ===================================================
% Fill-in the transformation equations for x(xi,eta) and y(xi,eta)
% for i = 1 : N_xi
%    for j = 1 : N_eta
%       xi  = (i-1)*d_xi;
%       eta = (j-1)*d_eta;
%       x(i,j)= ((cc-bb)*eta+bb)*xi;
%       y(i,j)= eta*hh;
%    end
% end


%%=== Solution output =============================================================
% Uncomment the following lines to make the plots
%
% %================================================================================
% %-- Plot the mesh
% figure,
% mesh(x,y,0*x,0*y);
% view([0,0,1]);
% axis equal
%
% %================================================================================
% %-- Plot the solution
% figure,
% [ccc,fff]=contour(x,y,abs(Solution),(0.02:0.07:max(max(abs(Solution))))');
% clabel(ccc,fff);
% hold on
% patch([0,bb,cc,0],[0,0,hh,hh],-ones(1,4),0,'facecolor',[0.8,0.8,0.8]);
% axis equal
