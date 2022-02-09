function [b] = interpolationFunctionDiff(phi,phimin,phimax)
%Set an derivertive of the interpolation function
b = max(3*phimin^2,min(3*phimax^2,3*phi^2));
% b = 3*phi^2;
% b=(phimin+(1-phimin))*3*phi^2;
end