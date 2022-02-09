function phaseField_inverter(h,nelx,nely,e,nu,lv,lag,kappa,eta,wcoef,PDECycle,timeStep,u_target_out)
ne = nelx*nely; %Number of elements
[~,nodes] = mesh(nelx,nely,h); %Set nodes position and nodes constructing each element.
phi = initialShape(nelx,nely); %Set initial value of phase field function phi.
phimin=0.01;
phimax=0.99;
%Start iteration
%while sum(vol)>=v0
for l = 1:10
    phi =  max(phimin,min(phimax,phi)); %Clip phi the range between 0 to 1. This is a charm against numerical errors.
    colormap(gray);
    or_im = flipud(reshape(-phi,nelx,nely)');
    fl_im = reshape(-phi,nelx,nely)';
    ful_im = imagesc(vertcat(fl_im,or_im),[-1.0 0]);
    axis equal; axis tight; axis off;pause(1e-6); %Above 3 lines are plotting phi.
    u = FEM(nelx,nely,nodes,e,nu,phi,phimin,phimax,lv); %Calculate displacement u by FEM.
    com = compliance(u,nelx,nely,u_target_out); %Calculate compliance
    dcom = complianceDiff(nelx,nely,nodes,e,nu,u,phi,phimin,phimax,u_target_out); %Calculate sensitivity of compliance
    vol = h^2 * phi;%Calculate volume
    dvol = h^2 * ones(ne,1);%Calculate sensitivity of volume
    objf = com + lag*sum(vol);%Calculate objective function
    disp([' It.: ' sprintf('%4i',l) ' Obj.: ' sprintf('%10.4f',objf) ' Vol.: ' sprintf('%6.3f',sum(vol)) ' Lag.: ' sprintf('%6.3f',lag)])
    dobjf = dcom + lag*dvol; %Calculate sensitivity
    avedobjf=sum(abs(dobjf))/length(dobjf);
    dobjf = dobjf/avedobjf; %Normalizing sensitivity value by its average
    G = eta*dobjf;
    W = wcoef;%Coefficient of function w(phi).
    
    for j = 1:PDECycle %Solve Allen-Cahn equation
        mphi =-((1 - 2 * phi) .* 2 * W + 30 * phi .* (1 - phi) .* G); %Calculate a part of reaction term, which is used in the semi-implicit method.
        phi = semiImplicit(phi,mphi,nelx,nely,h,timeStep,kappa);
    end
end

function [xyz,nodes] = mesh(nelx,nely,h)
% Generate an analysis domain and mesh. Nodes and Elements numbers are set as follows.
%
%11-12-13-14-15 %Node number
% | 5| 6| 7| 8| %Element number
% 6- 7- 8- 9-10 %Node number
% | 1| 2| 3| 4| %Element number
% 1- 2- 3- 4- 5 %Node number
%

xyz = zeros((nelx+1)*(nely+1),2); %Coordinates of each node
nodes = zeros(nelx*nely,4); %Index of nodes constructing each element

for i = 1:nely+1
    for j = 1:nelx+1
        xyz((i-1)*(nelx+1)+j,:) = [(j-1)*h (i-1)*h];
    end
end

for i = 1:nely
    for j = 1:nelx
        nodes((i-1)*nelx+j,:) = [i*(nelx+1)+j i*(nelx+1)+j+1 (i-1)*(nelx+1)+j+1 (i-1)*(nelx+1)+j];
    end
end

function [u] = FEM(nelx,nely,nodes,e,nu,phi,phimin,phimax,lv)
%Function of FEM analysis
ne = nelx*nely; %Number of elements
din = 2*(nelx+1)*(nely+1)-1;
dout = 2*(nelx+1)*nely+1;
ndofs = 2*(nelx+1)*(nely+1);
no_fixed_nodes = (0.1 / (1/nely))+1;
gk = sparse(ndofs, ndofs); %Global stiffness matrix
f = sparse(ndofs, 1); %Global load vector
f(din) = lv; %Substituting the local load vector into the global one
gk = globalStiffnessMatrix(ne,nodes,e,nu,gk,phi,phimin,phimax); %Construct global stiffness matrix
gk(din,din) = gk(din,din) + 0.01;
gk(dout,dout) = gk(dout,dout) + 0.01;
fixed_x_dofs = [2*(nelx+1)-1 :2*(nelx+1): 2*(nelx+1)*no_fixed_nodes-1];
fixed_y_dofs1 = [2*(nelx+1) :2*(nelx+1): 2*(nelx+1)*no_fixed_nodes];
fixed_y_dofs2 = [2*(nelx+1)*nely+2 :2: 2*(nelx+1)*(nely+1)];
fixedDofs = horzcat(fixed_x_dofs, fixed_y_dofs1, fixed_y_dofs2);
allDofs = [1:ndofs];
freeDofs = setdiff(allDofs,fixedDofs);
u(freeDofs,:) = gk(freeDofs,freeDofs) \ f(freeDofs,:); %Solve linear system
u(fixedDofs,:) = 0;


function [ke] = elementStiffnessMatrix(e,nu)

%Calculate element stiffness matrix. I reffered sigmund's 99 lines matrab
%topology optimization code.
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
ke = e/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

function [gk] = globalStiffnessMatrix(ne,nodes,e,nu,gk,phi,phimin,phimax)
%Construct global stiffness matrix
ke = elementStiffnessMatrix(e,nu);
for i = 1:ne
    nd = nodes(i,:);
    edofs = [nd(1)*2-1;nd(1)*2;nd(2)*2-1;nd(2)*2;nd(3)*2-1;nd(3)*2;nd(4)*2-1;nd(4)*2];
    gk(edofs,edofs) = gk(edofs,edofs) + interpolationFunction(phi(i),phimin,phimax)*ke;
end

function [com] = compliance(u,nelx,nely,u_target_out)
dout = 2*(nelx+1)*nely+1;
alpha = 3;
%Create target displacement vector
com = abs(u(dout)-u_target_out)^alpha;

function [dcom] = complianceDiff(nelx,nely,nodes,e,nu,u,phi,phimin,phimax,u_target_out)
%Calculate sensitivities of compliance
ne = nelx*nely; %Number of elements
din = 2*(nelx+1)*(nely+1)-1;
dout = 2*(nelx+1)*nely+1;
ndofs = 2*(nelx+1)*(nely+1);
no_fixed_nodes = (0.1 / (1/nely))+1;
dcom = zeros(ne,1);
u_target = u_target_out;
alpha = 3;

ke = elementStiffnessMatrix(e,nu);
gk = sparse(ndofs, ndofs);
gk = globalStiffnessMatrix(ne,nodes,e,nu,gk,phi,phimin,phimax);
gk(din,din) = gk(din,din) + 0.01;
gk(dout,dout) = gk(dout,dout) + 0.01;
f_load = sparse(ndofs, 1);
du = u(dout,1) - u_target;
C0 = (abs(du)^alpha)^(1/(1/alpha));
f_load(dout,1) = C0* abs(du)^(alpha-2);

allDofs = [1:ndofs];
fixed_x_dofs = [2*(nelx+1)-1 :2*(nelx+1): 2*(nelx+1)*no_fixed_nodes-1];
fixed_y_dofs1 = [2*(nelx+1) :2*(nelx+1): 2*(nelx+1)*no_fixed_nodes];
fixed_y_dofs2 = [2*(nelx+1)*nely+2 :2: 2*(nelx+1)*(nely+1)];
fixedDofs = horzcat(fixed_x_dofs, fixed_y_dofs1, fixed_y_dofs2);
freeDofs = setdiff(allDofs,fixedDofs);
p(freeDofs,:) = gk(freeDofs,freeDofs) \ (f_load(freeDofs,:)); %Solve linear system
p(fixedDofs,:) = 0;

for i=1:ne
    nd = nodes(i,:);
    edofs = [nd(1)*2-1;nd(1)*2;nd(2)*2-1;nd(2)*2;nd(3)*2-1;nd(3)*2;nd(4)*2-1;nd(4)*2];
    dcom(i) = interpolationFunctionDiff(phi(i),phimin,phimax)*u(edofs,:)'*ke*p(edofs,:);
end


function [b] = interpolationFunction(phi,phimin,phimax)
%Set an interpolation function of diffuse interface domain
b = max(phimin,min(phimax,phi^3));

function [b] = interpolationFunctionDiff(phi,phimin,phimax)
%Set an derivertive of the interpolation function
b = max(3*phimin^2,min(3*phimax^2,3*phi^2));

function [phi] = semiImplicit(phi,mphi,nelx,nely,h,dt,eps)

%Semi-implicit method for Allen-Cahn equation.
%Please refer to J. A. Warren, R. Kobayashi, A. E. Lobkovsky, and W. C. Carter, Acta Mater. 51, 6035 (2003)
%Warning! Since the code is not optimized for matrab, the calcuration is quite slow despite nearly explicit method!
% dx = h ; dy = h ;
newphi = zeros(nely,nelx);
phi = reshape(phi,nelx,nely)'; %Reshape phase field function for easy calcuration of finite difference method.
mphi = reshape(mphi,nelx,nely)';

for i = 1:nelx
    for j = 1:nely
        if mphi(j,i) <= 0
            tmp = mphi(j,i)*(1-phi(j,i));
            if 1<i && i<nelx
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            elseif i == 1
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            else
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            end
        else
            tmp = mphi(j,i)*phi(j,i);
            if 1<i && i<nelx
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i))  + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i))  + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i))  + dt*tmp)/(1 + dt*tmp);
                end
            elseif i == 1
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - phi(j,i)) + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(phi(j,i+1) - phi(j,i) - phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                end
            else
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)) + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (eps*dt/h^2)*(- phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                end
            end
        end
    end
end
phi = reshape(newphi,nely,nelx)';
phi = reshape(phi,nelx*nely,1);

function [phi] = initialShape(nelx,nely)
%Make an initial shape with several holes.
phi = ones(nelx,nely);
r = nely/8.0;%Radius of holes
rc = [nelx/6.0 0; nelx/2.0 0; nelx*5.0/6.0 0;
    0 nely/4.0; nelx/3.0 nely/4.0; nelx*2.0/3.0 nely/4.0; nelx nely/4.0;
    nelx/6.0 nely/2.0; nelx/2.0 nely/2.0; nelx*5.0/6.0 nely/2.0];

%Center of holes
for l = 1:length(rc)
    center = rc(l,:);
    for i = 1:nelx
        for j = 1:nely
            if norm(center - [i j])<r
                phi(i,j) = 0;
                phi(i,nely-j+1) = 0;
            end
        end
    end
end
phi = reshape(phi,nelx*nely,1);