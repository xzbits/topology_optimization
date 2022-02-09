function [b] = interpolationFunction(phi,phimin,phimax)
%Set an interpolation function of diffuse interface domain
b = max(phimin,min(phimax,phi^3));
% b=(phimin+(1-phimin))*phi^3;
end