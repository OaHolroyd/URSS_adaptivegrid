function V = periodicInterp(Xold,X,V,L,interpmethod)
% slight adaptation of MATLAB's standard 'interp1' which better accounts for
% periodic boundaries

wrap = 5; % number of points to pad the ends with

% add duplicates to account to periodicity
Xold = [Xold(end-wrap:end)-L;Xold;Xold(1:wrap)+L];
V = [V(end-wrap:end);V;V(1:wrap)];

V = interp1(Xold,V,X,interpmethod);
end
