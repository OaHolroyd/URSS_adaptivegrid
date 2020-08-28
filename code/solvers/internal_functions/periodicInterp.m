function V = periodicInterp(Xold,X,V,L,interpmethod)
%PERIODICINTERP 'interp1' with periodic boundaries
% slight adaptation of MATLAB's standard 'interp1' which better accounts for
% periodic boundaries
%
%   INPUTS: Xold - origional grid (for V)
%              X - grid V is to be interpolated to
%              V - vector to be interpolated
%              L - length of domain. Should be greater than X(end)
%   interpmethod - method used to interpolate from Xold to X
%
%
%   OUTPUTS: V - V interpolated to X

wrap = 5; % number of points to pad the ends with

% add duplicates to account to periodicity
Xold = [Xold(end-wrap:end)-L;Xold;Xold(1:wrap)+L];
V = [V(end-wrap:end);V;V(1:wrap)];

V = interp1(Xold,V,X,interpmethod);
end
