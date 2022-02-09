%Matlab code of the phase field method for shape optimization by Akihiro
%Takezawa 2016.09.16
clc
clear all

e = 1.0; %Young's modulus
nu = 0.3; %Poisson's ratio
lx = 1.0; %Length of x
ly = 1.0; %Length of y
nelx = 100; %Discritization number of x
nely = 100; %Discritization number of y

% This code only assumes only square mesh.
h = lx / nelx; %Element size
lv = 1; %Local load vector
lp = (nelx+1)*1; %load point nodes number (center of the right side)
lag = 0.5; %Lagrange multiplier for volume constraint
u_target_out = 5.6;

% Compliant mechanism
a = 1; %Parameter of double well potential representing right side height. It should not be changed.
thickness = 9*h; % Thickness of phase field interface with 0.05<phi<0.95. Can be changed.
ep = thickness*a/6; % Square root of diffusion coefficient.

kappa = 1e-5; % Diffusion coefficient.
eta = 1; %Coefficient of sensitivity, 0.5-1.0 works well.
wcoef = 0.25; %Parameter of double well potential representing center bump height, it should not be changed.

PDECycle = 1; %Number of the time updating in Allen-Cahn equation per one iteration.
%In some sensitive problems (such as compliant mechanism problem or vibration problem), it should be small.
timeStep = 0.9 * h^2/(4*kappa); %It must satisfy CFL condition.

phaseField_auxetic(h,nelx,nely,e,nu,lv,lag,kappa,eta,wcoef,PDECycle,timeStep,u_target_out);