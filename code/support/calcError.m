function err = calcError(X,U,Xsol,Usol,L)
%CALCERROR Calculates the 2-norm error between U and Usol

V = periodicInterp(Xsol,X,Usol,L,'pchip');

err = sqrt(trapz(X,(abs(U-V).^2)));
end

