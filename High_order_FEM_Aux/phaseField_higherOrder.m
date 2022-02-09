function phaseField_higherOrder(h,L,elemType,nelx,nely,C,lv,lp,lag,kappa,eta,wcoef,PDECycle,timeStep,u_target_out,n_target_disp,uln,lrn)
ne = nelx*nely; %Number of elements
phi = initialShape(nelx,nely); %Set initial value of phase field function phi.
phimin=0.01;
phimax=0.99;
%Start iteration
for i=1:40
    phi =  max(phimin,min(phimax,phi)); %Clip phi the range between 0 to 1. This is a charm against numerical errors.
    colormap(gray);
    
    %Visualize the optimal shape
    or_im = flipud(reshape(-phi,nelx,nely)');
    fl_im = reshape(-phi,nelx,nely)';
    h1_im = vertcat(or_im,fl_im);
    h2_im = flipdim(h1_im, 2);
    ful_im = imagesc(horzcat(h2_im,h1_im),[-1.0 0]);
    
    axis equal; axis tight; axis off;pause(1e-6); %Above 3 lines are plotting phi.
    [u,nodes,~] = FEM(L,elemType,nelx,nely,C,lv,lp,phi,phimin,phimax); %Calculate displacement u by FEM.
    numnode=size(nodes,1);
    com = compliance(u,numnode,u_target_out,uln); %Calculate compliance
    dcom = complianceDiff(L,elemType,nelx,nely,numnode,C,phi,phimin,phimax,u,u_target_out,n_target_disp,uln,lrn); %Calculate sensitivity of compliance
    vol = h^2 * phi;%Calculate volume
    dvol = h^2 * ones(ne,1);%Calculate sensitivity of volume
    objf = com + lag*(sum(vol));%Calculate objective function
    disp([' It.: ' sprintf('%4i',i) ' Obj.: ' sprintf('%10.4f',objf) ' Vol.: ' sprintf('%6.3f',sum(vol)) ' u_in: ' sprintf('%10.4f',u(lrn)) ' u_out: ' sprintf('%10.4f',u(uln+numnode)) ' poisson: ' sprintf('%10.4f',-u(uln+numnode)/u(lrn)) ' lv: ' sprintf('%10.4f',lv(1))])
    dobjf = dcom + lag*(dvol); %Calculate sensitivity
    avedobjf=sum(abs(dobjf))/length(dobjf);
    dobjf = dobjf/avedobjf; %Normalizing sensitivity value by its average
    G = eta*dobjf;
    W = wcoef;%Coefficient of function w(phi).
    for j = 1:PDECycle %Solve Allen-Cahn equation
        mphi =-((1 - 2 * phi) .* 2 * W + 30 * phi .* (1 - phi) .* G); %Calculate a part of reaction term, which is used in the semi-implicit method.
        phi = semiImplicit(phi,mphi,nelx,nely,h,timeStep,kappa);
    end
   
end %end volume loop

function [com] = compliance(u,numnode,u_target_out,uln)
dout = uln + numnode;
alpha = 3;
%Create target displacement vector
com = abs(u(dout)-u_target_out)^alpha;

function [dcom] = complianceDiff(L,elemType,nelx,nely,numnode,C,phi,phimin,phimax,u,u_target_out,n_target_disp,uln,lrn)
%Calculate sensitivities of compliance
numelem=nelx*nely;
din = lrn;
dout = uln + numnode;

alpha = 3;
du = u(dout,1) - u_target_out;
C0 = (abs(du)^alpha)^(1/(1/alpha));
f_load = C0* abs(du)^(alpha-2);

[p,nodes,elements] = FEM(L,elemType,nelx,nely,C,f_load,n_target_disp,phi,phimin,phimax);

switch elemType % define quadrature rule
case 'Q9'
[W,Q]=quadrature( 4, 'GAUSS', 2 ); % 4x4 Gaussian quadrature
case 'Q4'
[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
end

dcom=zeros(numelem,1);
for el=1:numelem % start of element loop
    sctr=elements(el,:); % element scatter vector
    sctrB=[ sctr sctr+numnode ]; % vector that scatters a B matrix
    nn=length(sctr);
    K=sparse(length(sctrB),length(sctrB));
    for q=1:size(W,1) % quadrature loop
        pt=Q(q,:); % quadrature point
        wt=W(q); % quadrature weight
        [~,dNdxi]=shape_func(elemType,pt); % element shape functions
        J0=nodes(sctr,:)'*dNdxi; % element Jacobian matrix
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K=K+B'*C*B*wt*det(J0);       
    end % of quadrature loop
    dcom(el) = u(sctrB)'*interpolationFunctionDiff(phi(el),phimin,phimax)*K*p(sctrB);
end % of element loop

function [phi] = semiImplicit(phi,mphi,nelx,nely,h,dt,kappa)

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
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            elseif i == 1
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            else
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)))/(1 - dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)))/(1 - dt*tmp);
                end
            end
        else
            tmp = mphi(j,i)*phi(j,i);
            if 1<i && i<nelx
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i))  + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i))  + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - 2*phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i))  + dt*tmp)/(1 + dt*tmp);
                end
            elseif i == 1
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) + phi(j+1,i) - phi(j,i)) + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(phi(j,i+1) - phi(j,i) - phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                end
            else
                if 1<j && j<nely
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
                elseif j == 1
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) + phi(j+1,i) - phi(j,i)) + dt*tmp)/(1 + dt*tmp);
                else
                    newphi(j,i) = (phi(j,i) + (kappa*dt/h^2)*(- phi(j,i) + phi(j,i-1) - phi(j,i) + phi(j-1,i)) + dt*tmp)/(1 + dt*tmp);
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
r = nely/12.0;%Radius of holes
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