function [U,node,element] = FEM(L,elemType,numx,numy,C,lv,lp,phi,phimin,phimax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE MESH
%
% node - is a matrix of the node coordinates, i.e. node(I,j) -> x_Ij
% element - is a matrix of element connectivities, i.e. the connectivity
% of element e is given by > element(e,:) -> [n1 n2 n3 ...];
%
% rightEdge - a element connectivity matrix for the right edge
% leftEdge
%
% These connectivity matricies refer to the node numbers defined in the
% coordinate matrix node.
switch elemType
case 'Q4' % here we generate the mesh of Q4 elements
nnx=numx+1;
nny=numy+1;
node=square_node_array([0 0],[L 0],[L L],[0 L],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
case 'Q9' % here we generate a mehs of Q9 elements
nnx=2*numx+1;
nny=2*numy+1;
node=square_node_array([0 0],[L 0],[L L],[0 L],nnx,nny);
inc_u=2;
inc_v=2*nnx;
node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];
element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
end
% DEFINE BOUNDARIES
uln=nnx*(nny-1)+1; % upper left node number
urn=nnx*nny; % upper right node number
lrn=nnx; % lower right node number
lln=1; % lower left node number
% GET NODES ON DISPLACEMENT BOUNDARY
% Here we get the nodes on the essential boundaries
fixedNodeX=[lln:nnx:uln]'; % a vector of the node numbers which are fixed in
% the x direction
fixedNodeY=[lln:lrn]'; % a vector of node numbers which are fixed in
% the y-direction

uFixed=zeros(size(fixedNodeX)); % a vector of the x-displacement for fixedNodeX
vFixed=zeros(size(fixedNodeY)); % and the y-displacements for fixedNodeY

numnode=size(node,1); % number of nodes
numelem=size(element,1); % number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES
%
% Here we define the system data structures
% U - is vector of the nodal displacements it is of length 2*numnode. The
% displacements in the x-direction are in the top half of U and the
% y-displacements are in the lower half of U, for example the displacement
% in the y-direction for node number I is at U(I+numnode)
% f - is the nodal force vector. It's structure is the same as U,
% i.e. f(I+numnode) is the force in the y direction at node I
% K - is the global stiffness matrix and is structured the same as with U and f
% so that K_IiJj is at K(I+(i-1)*numnode,J+(j-1)*numnode)
f=zeros(2*numnode,1); % external load vector
K=sparse(2*numnode,2*numnode); % stiffness matrix
% ******************************************************************************
% *** P R O C E S S I N G ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the force at the nodes on the top and bottom edges to zero (traction free)
% TOP EDGE
% topEdgeNodes = find(node(:,2)==L); % finds nodes on the top edge
% f(topEdgeNodes)=0;
% f(topEdgeNodes+numnode)=0;
% % RIGHT EDGE
% rightEdgeNodes = find(node(:,1)==L); % finds nodes on the bottom edge
% f(rightEdgeNodes)=0;
% f(rightEdgeNodes+numnode)=0;
% FORCE ON NODE
switch lp
case 'load'
f(lrn)=lv;
case 'target disp'
f(uln+numnode)=lv;
end
%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch elemType % define quadrature rule
case 'Q9'
[W,Q]=quadrature( 4, 'GAUSS', 2 ); % 4x4 Gaussian quadrature
case 'Q4'
[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
end
for e=1:numelem % start of element loop
    Cphi=interpolationFunction(phi(e),phimin,phimax)*C;
    sctr=element(e,:); % element scatter vector
    sctrB=[ sctr sctr+numnode ]; % vector that scatters a B matrix
%     sctrB=[ sctr(1)*2-1 sctr(1)*2 sctr(2)*2-1 sctr(2)*2 sctr(3)*2-1 sctr(3)*2 sctr(4)*2-1 sctr(4)*2 ]; % vector that scatters a B matrix
%     sctrB=[sctr(4)*2-1 sctr(4)*2 sctr(3)*2-1 sctr(3)*2 sctr(2)*2-1 sctr(2)*2 sctr(1)*2-1 sctr(1)*2];
    nn=length(sctr);
    kel=zeros(2*nn,2*nn);
    for q=1:size(W,1) % quadrature loop
        pt=Q(q,:); % quadrature point
        wt=W(q); % quadrature weight
        [~,dNdxi]=shape_func(elemType,pt); % element shape functions
        J0=node(sctr,:)'*dNdxi; % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE B MATRIX
        % _ _
        %       | N_1,x     N_2,x   ... 0       0 ...       |
        % B =   | 0         0       ... N_1,y   N_2,y ...   |
        %       | N_1,y     N_2,y   ... N_1,x   N_2,x ...   |
        % - -
        B=zeros(3,2*nn);
        B(1,1:nn) = dNdx(:,1)';
        B(2,nn+1:2*nn) = dNdx(:,2)';
        B(3,1:nn) = dNdx(:,2)';
        B(3,nn+1:2*nn) = dNdx(:,1)';

%         B=zeros(3,2*nn);
%         B(1,1:2:2*nn) = dNdx(:,1)';
%         B(2,2:2:2*nn) = dNdx(:,2)';
%         B(3,1:2:2*nn) = dNdx(:,2)';
%         B(3,2:2:2*nn) = dNdx(:,1)';


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT    
%         K(sctrB,sctrB)=K(sctrB,sctrB)+B'*Cphi*B*wt*det(J0);
        kel=kel+B'*C*B*wt*det(J0);
    end % of quadrature loop
    K(sctrB,sctrB)=K(sctrB,sctrB)+interpolationFunction(phi(e),phimin,phimax)*kel;
end % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%
% APPLY ESSENTIAL BOUNDARY CONDITIONS
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix

udofs=fixedNodeX; % global indecies of the fixed x displacements
vdofs=fixedNodeY+numnode; % global indecies of the fixed y displacements

% udofs=2*fixedNodeX-1; % global indecies of the fixed x displacements
% vdofs=2*fixedNodeY; % global indecies of the fixed y displacements
din = lrn;
dout = uln + numnode;
f=f-K(:,udofs)*uFixed; % modify the force vector
f=f-K(:,vdofs)*vFixed;
f(udofs)=bcwt*uFixed;
f(vdofs)=bcwt*vFixed;
K(udofs,:)=0; % zero out the rows and columns of the K matrix
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(din,din) = K(din,din) + 0.01;
K(dout,dout) = K(dout,dout) + 0.01;
% SOLVE SYSTEM
U=K\f; % nodal displacement vector

%******************************************************************************
%*** P O S T - P R O C E S S I N G ***
%******************************************************************************

% a vector of indicies that quickly address the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal x-displacements
% xs=1:numnode; % x portion of u and v vectors
% ys=(numnode+1):2*numnode; % y portion of u and v vectors
% 
% dispNorm=L/max(sqrt(U(xs).^2+U(ys).^2));
% scaleFact=0.1*dispNorm;
% 
% figure
% plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,U(xs));
% % hold on
% % plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
% % plot_mesh(node,element,elemType,'w--');
% colorbar
% title('DEFORMED DISPLACEMENT IN X-DIRECTION')
% 
% figure
% plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,U(ys));
% colorbar
% title('DEFORMED DISPLACEMENT IN Y-DIRECTION')
end