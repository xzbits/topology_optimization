%Matlab code of the phase field method for shape optimization by Akihiro
%Takezawa 2016.09.16
clc
clear all

e = 1; %Young's modulus
nu = 0.3; %Poisson's ratio

nelx = 100; %Discritization number of x
nely = 100; %Discritization number of y

L = 1.0; 
% This code only assumes only square mesh.
h = L / nelx;

elemType = 'Q4';
switch elemType
case 'Q4' % here we generate the mesh of Q4 elements
nnx=nelx+1;
nny=nely+1;
case 'Q9' % here we generate a mehs of Q9 elements
nnx=2*nelx+1;
nny=2*nely+1;
end

% DEFINE BOUNDARIES
uln=nnx*(nny-1)+1; % upper left node number
urn=nnx*nny; % upper right node number
lrn=nnx; % lower right node number
lln=1; % lower left node number

load = 1; %Local load vector
u_target_out = 7.0;

% STRESS ASSUMPTION
stressState='PLANE_STRESS'; % set to either 'PLANE_STRAIN' or "PLANE_STRESS'
% COMPUTE ELASTICITY MATRIX
if ( strcmp(stressState,'PLANE_STRESS') ) % Plane Stress case
C=e/(1-nu^2)*[ 1 nu 0;
nu 1 0;
0 0 (1-nu)/2 ];
else % Plane Strain case
C=e/(1+nu)/(1-2*nu)*[ 1-nu nu 0;
nu 1-nu 0;
0 0 1/2-nu ];
end

% Compliant mechanism
lag = 0.5; %Lagrange multiplier for volume constraint
a = 1; %Parameter of double well potential representing right side height. It should not be changed.
thickness = 9*h; % Thickness of phase field interface with 0.05<phi<0.95. Can be changed.
ep = thickness*a/6; % Square root of diffusion coefficient.

kappa = ep^2; % Diffusion coefficient.
eta = 1; %Coefficient of sensitivity, 0.5-1.0 works well.
wcoef = 0.25; %Parameter of double well potential representing center bump height, it should not be changed.

PDECycle = 3; %Number of the time updating in Allen-Cahn equation per one iteration.
%In some sensitive problems (such as compliant mechanism problem or vibration problem), it should be small.
timeStep = 0.3 * h^2/(4*kappa); %It must satisfy CFL condition.

phaseField_higherOrder(h,L,elemType,nelx,nely,C,load,'load',lag,kappa,eta,wcoef,PDECycle,timeStep,u_target_out,'target disp',uln,lrn);