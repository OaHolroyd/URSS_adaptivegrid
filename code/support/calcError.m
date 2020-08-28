function err = calcError(X,U,Xsol,Usol,L)
%CALCERROR Calculates the 2-norm error between U and Usol
%
%   INPUTS:  X - grid for U
%            U - vector to be compared
%         Xsol - grid for Usol
%         Usol - vector to be compared to
%            L - length of domain
%
%   OUTPUTS: err - error

V = periodicInterp(Xsol,X,Usol,L,'pchip');

err = sqrt(trapz(X,(abs(U-V).^2)));
end

