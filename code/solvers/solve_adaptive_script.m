%% IMPLICIT NEWT ADAPTIVE BDF2 EQUIMESH
% Solves the Benney IVP on an adaptive grid. Uses BDF1 for the first timestep
% and the does the rest using BDF2. Uses 2nd order spatial derivative operators.
% Adaptive grid with padding as in Kautsky & Nichols 1980
%  INPUTS: bparams - Benney PARAMeterS. Path to a .mat file containing the
%                    solver-independent parameters:
%                    g -
%          mparams - Model PARAMeterS. Path to a .mat file containing the
%                    solver-dependent parameters.
%
% OUTPUTS: Xout - non-uniform grids at Tout
%          Uout - interface height at Xout/Tout
%          Fout - forcing term at Xout/Tout
%          Wout - weighting function at Xout/Tout
%          Pout - padded weight at Xout/Tout
%  elasped_time - CPU time spent on solving (doesn't include setup)
%             N - number of gridpoints used

%% Paramters
% load params from file
load(bparams,'g','f','L','Tmax','b','R','C');
load(mparams,'dxinf','dxsup','K','d','Tout','interpmethod','filename');

MAXITER = 50; % maximum number of iterations for Newton/Broyden solver
ITERACC = 5*10^-9; % threshhold for the solver to stop

GRIDUPDATE = 0; % how frequently to update the grid
c = 1; % dt/dx

cotb = cot(b); % precompute cot(b);
Ntotal = 0;
Ncount = 0;


%% Storage
Tout = sort(Tout);
if Tout(1)~=0
    Tout = [0;Tout];
end
Mout = length(Tout);

Uout = cell(Mout,1); % for interface height
Fout = cell(Mout,1); % for controls
Xout = cell(Mout,1); % for the grids
Wout = cell(Mout,1); % for the weights
Pout = cell(Mout,1); % for the padded weights

t = 0; % current time
dt_1 = c*dxinf; % previous timestep
dt0 = c*dxinf; % current timestep
jmax = 1; % highest output index reached (in case of early finish)\


%% Uniform precomputing
Nuni = 512;
Xuni = linspace(0,L,Nuni+1)';
Xuni = Xuni(1:end-1);
dxuni = Xuni(2) - Xuni(1);

% uniform grid derivatives
Dxuni = derMat(1,Nuni,dxuni);
Dxxuni = derMat(2,Nuni,dxuni);

U0uni = g(Xuni); % change to vector at some point probably
F0uni = f(Xuni,U0uni,0); % change to vector at some point probably


%% Non-uniform initial settup
W = calcWeight(U0uni,F0uni,Dxuni,Dxxuni);
X = equiMesh(Xuni,W,L,K,d,dxinf,dxsup);
N = length(X);
I = speye(N,N);
Ntotal = Ntotal + N; Ncount = Ncount + 1;

% interpolate to new grid
U = periodicInterp(Xuni,X,U0uni,L,interpmethod);
U0 = U;
U_1 = U;
F = periodicInterp(Xuni,X,F0uni,L,interpmethod);

% repeat
D1 = derMatNonUniform(1,2,X,L);
D2 = derMatNonUniform(2,3,X,L);
D3 = derMatNonUniform(3,2,X,L);
D4 = derMatNonUniform(4,3,X,L);

Xold = X;
W = calcWeight(U,F,D1,D2);
[X,P,dxmin] = equiMesh(X,W,L,K,d,dxinf,dxsup);
N = length(X);
I = speye(N,N);
dt = c*dxmin; % set dt from dxmin

U_1 = periodicInterp(Xold,X,U_1,L,interpmethod);
U0 = periodicInterp(Xold,X,U0,L,interpmethod);
U = periodicInterp(Xold,X,U,L,interpmethod);
F = periodicInterp(Xold,X,F,L,interpmethod);

% store initial conditions
if Tout(1)==0
    Xout{1} = X;
    Uout{1} = U;
    Fout{1} = F;
    Wout{1} = W;
    Pout{1} = P;
end

%% First step (BDF1)
tic; % ignore setup when timing

% check if an output is required. If so, adjust dt accordingly
output = false;
if (t + dt)>Tout(jmax+1)
    output = true;
    dt0 = Tout(jmax+1)-t;
elseif (t + dt)==Tout(jmax+1)
    output = true;
    dt0 = dt;
else
    dt0 = dt;
end

t = t+dt0;

F = f(X,U,t); % update control vector

% update grid/interp to new grid: already done so skip this on the first step
Ntotal = Ntotal + N; Ncount = Ncount + 1;

D1 = derMatNonUniform(1,2,X,L);
D2 = derMatNonUniform(2,3,X,L);
D3 = derMatNonUniform(3,2,X,L);
D4 = derMatNonUniform(4,3,X,L); % derivatives on non-uniform grid

% Newton/Broyden solver
G = calcG(U,F,D1,D2,D3,D4,R,C,cotb);
fU = - dt0 * G;

Jg = calcJg(U,F,I,D1,D2,D3,D4,R,C,cotb); % Jacobian for G
J = I - dt0 * Jg;

% iterate solver (for at most MAXITER iterations)
for k = 1:MAXITER
    delU = J\-fU;
    U = U + delU;
    G = calcG(U,F,D1,D2,D3,D4,R,C,cotb);
    delfU = (U - U0 - dt0 * G) - fU;
    fU = fU + delfU;

    % stop iterating if close to a solution
    if norm(fU) < ITERACC
        break
    end

    % good Broyden method to update J
    J = J + (delfU - J*delU)/sum(delU.*delU) * delU';
end

% abort if something has gone wrong
if sum(isnan(U)) > 0
    elapsed_time = toc;

    %remove unused storage
    Tout = Tout(1:jmax);
    Uout = Uout(1:jmax);
    Fout = Fout(1:jmax);
    Xout = Xout(1:jmax);
    Wout = Wout(1:jmax);
    Pout = Pout(1:jmax);

    warning('NaN value occured.\n%s\n%s\n%s',...
        ['Aborted at t = ',num2str(t)],...
        [num2str(t*100/Tmax,3),'% completed'],...
        ['jmax = ',num2str(jmax)]);
    %save(['saved_data/',filename],'Tout','Xout','Uout','Fout');
    return
end

% output
if output
    jmax = jmax + 1;
    Uout{jmax} = U;
    Fout{jmax} = F;
    Xout{jmax} = X;
    Wout{jmax-1} = W;
    Pout{jmax-1} = P;
end

% store previous iterations
U_1 = U0;
U0 = U;
dt_1 = dt0;


%% Remaining steps (BDF2)
while t<=Tmax
    F = f(X,U,t); % preliminary control vector update

    % update the grid if required
    if isnan(rem(t,GRIDUPDATE)) || rem(t,GRIDUPDATE)<dt
        % update grid
        Xold = X;
        %W = 0.5*(abs(Dx*U) + abs(Dx*F) + abs(Dxx*U) + abs(Dxx*F));
        W = calcWeight(U,F,D1,D2);
        [X,P,dxmin] = equiMesh(X,W,L,K,d,dxinf,dxsup);
        N = length(X);
        I = speye(N,N);

        % interp to new grid
        U_1 = periodicInterp(Xold,X,U_1,L,interpmethod);
        U0 = periodicInterp(Xold,X,U0,L,interpmethod);
        U = U0;
        F = periodicInterp(Xold,X,F,L,interpmethod);

        D1 = derMatNonUniform(1,2,X,L);
        D2 = derMatNonUniform(2,3,X,L);
        D3 = derMatNonUniform(3,2,X,L);
        D4 = derMatNonUniform(4,3,X,L); % derivatives on non-uniform grid

        % set timestep from dxmin
        dt = c*dxmin;
    end
    Ntotal = Ntotal + N; Ncount = Ncount + 1;

    % check if an output is required. If so, adjust dt accordingly
    output = false;
    if (t + dt)>Tout(jmax+1)
        output = true;
        dt0 = Tout(jmax+1)-t;
    elseif (t + dt)==Tout(jmax+1)
        output = true;
        dt0 = dt;
    else
        dt0 = dt;
    end

    t = t+dt0;


    % Newton/Broyden solver
    G = calcG(U,F,D1,D2,D3,D4,R,C,cotb);
    fU = (-dt0/dt_1/(dt0+dt_1))*U0 + (dt0/dt_1/(dt0+dt_1))*U_1 - G;

    Jg = calcJg(U,F,I,D1,D2,D3,D4,R,C,cotb); % Jacobian for G
    J = (2*dt0 + dt_1/dt0/(dt0+dt_1))*I - Jg;

    % iterate solver (for at most MAXITER iterations)
    for k = 1:MAXITER
        delU = J\-fU;
        U = U + delU;
        G = calcG(U,F,D1,D2,D3,D4,R,C,cotb);
        delfU = ((2*dt0 + dt_1)/dt0/(dt0+dt_1))*U + (-(dt0+dt_1)/dt0/dt_1)*U0 + (dt0/dt_1/(dt0+dt_1))*U_1 - G - fU;
        fU = fU + delfU;

        % stop iterating if close to a solution
        if norm(fU) < ITERACC
            break
        end

        % good Broyden method to update J
        J = J + (delfU - J*delU)/sum(delU.*delU) * delU';
    end

    % abort if something has gone wrong
    if sum(isnan(U)) > 0
        elapsed_time = toc;

        %remove unused storage
        Tout = Tout(1:jmax);
        Uout = Uout(1:jmax);
        Fout = Fout(1:jmax);
        Xout = Xout(1:jmax);
        Wout = Wout(1:jmax);
        Pout = Pout(1:jmax);

        warning('NaN value occured.\n%s\n%s\n%s',...
            ['Aborted at t = ',num2str(t)],...
            [num2str(t*100/Tmax,3),'% completed'],...
            ['jmax = ',num2str(jmax)]);
        %save(['saved_data/',filename],'Tout','Xout','Uout','Fout');
        return
    end

    % output
    if output
        jmax = jmax + 1;
        Uout{jmax} = U;
        Fout{jmax} = F;
        Xout{jmax} = X;
        Wout{jmax-1} = W;
        Pout{jmax-1} = P;
    end

    % store previous iterations
    U_1 = U0;
    U0 = U;
    dt_1 = dt0;

    if t>=Tout(end) || jmax>=Mout
        break
    end
end
elapsed_time = toc;

try
    W = calcWeight(U,F,D1,D2);
    [X,P,dxmin] = equiMesh(X,W,L,K,d,dxinf,dxsup);

    Wout{jmax} = W;
    Pout{jmax} = P;
catch
end

%% Function definitions
% computes the weighting vector used to set the grid
function W = calcWeight(U,F,D1,D2)
W = 0.5*(abs(D1*U) + 0.001*abs(D1*F) + abs(D2*U) + 0.001*abs(D2*F));

end

% computes G (where G(H,F) = Ht = F - Qx)
function G = calcG(U,F,D1,D2,D3,D4,R,C,cotb)
U2 = U.*U;
U3 = U2.*U;
U4 = U3.*U;
U5 = U4.*U;
U6 = U5.*U;

D1U = D1*U;
D2U = D2*U;
D3U = D3*U;
D4U = D4*U;
D1F = D1*F;

G = F - 2 * U2.*D1U...
    + 2*cotb * U2.*D1U.*D1U...
    + 2*cotb/3 * U3.*D2U...
    - 1/C * U2.*D1U.*D3U...
    - 1/3/C * U3.*D4U...
    - 16*R/5 * U5.*D1U.*D1U...
    - 8*R/15 * U6.*D2U...
    + 8*R/3 * U3.*D1U.*F...
    + 2*R/3 * U4.*D1F;
end

% computes the Jacobian for G for the Newton solver
function Jg = calcJg(U,F,I,D1,D2,D3,D4,R,C,cotb)
U2 = U.*U;
U3 = U2.*U;
U4 = U3.*U;
U5 = U4.*U;
U6 = U5.*U;

J2 = 2 * (I.*U);
J3 = 3 * (I.*U2);
J4 = 4 * (I.*U3);
J5 = 5 * (I.*U4);
J6 = 6 * (I.*U5);

D1U = D1*U;
D2U = D2*U;
D3U = D3*U;
D4U = D4*U;
D1F = D1*F;

Jg = -2 * ( J2.*D1U + D1.*U2 )...
     + 2*cotb * ( J2.*D1U.*D1U + 2*D1.*D1U.*U2 )...
     + 2*cotb/3 * ( J3.*D2U + D2.*U3 )...
     - 1/C * ( J2.*D1U.*D3U + D1.*D3U.*U2 + D3.*D1U.*U2 )...
     - 1/3/C * ( J3.*D4U + D4.*U3 )...
     - 16*R/5 * ( J5.*D1U.*D1U + 2*D1.*D1U.*U5 )...
     - 8*R/15 * ( J6.*D2U + D2.*U6 )...
     + 8*R/3 * ( J3.*D1U + D1.*U3 ) .* F...
     + 2*R/3 * ( J4.*D1F );
end

