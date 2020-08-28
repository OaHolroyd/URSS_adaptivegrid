%% SOLVE UNIFORM
% Solves the Benney IVP on a uniform grid. Uses BDF1 for the first timestep
% and the does the rest using BDF2. Uses 2nd order spatial derivative operators.
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
load(mparams,'dx','dt','Tout','filename');

N = floor(L/dx);

MAXITER = 10; % maximum number of iterations for Newton solver
ITERACC = 5 * 10^-9; % threshhold for the solver to stop

cotb = cot(b); % precompute cot(b);


%% Storage
% spatial grid
X = linspace(0,L,N+1)';
X = X(1:end-1);
dx = X(2) - X(1);

Tout = sort(Tout);
if Tout(1)~=0
    Tout = [0;Tout];
end
Mout = length(Tout);

% output storage
Uout = cell(Mout,1); % for interface height
Fout = cell(Mout,1); % for controls
Xout = cell(Mout,1); % for the grids

t = 0; % current time
dt_1 = dt; % previous timestep
dt0 = dt; % current timestep
jmax = 1; % highest output index reached (in case of early finish)


%% Uniform precomputing
% uniform grid derivatives
I = speye(N,N);
D1 = derMat(1,N,dx);
D2 = derMat(2,N,dx);
D3 = derMat(3,N,dx);
D4 = derMat(4,N,dx);

U = g(X); % change to vector at some point probably
U0 = U;
U_1 = U;
F = f(X,U,t); % change to vector at some point probably

% store initial conditions if required
if Tout(1)==0
    Xout{1} = X;
    Uout{1} = U;
    Fout{1} = F;
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

    % remove unused storage
    Tout = Tout(1:jmax);
    Uout = Uout(1:jmax);
    Fout = Fout(1:jmax);
    Xout = Xout(1:jmax);

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
end

% store previous iterations
U_1 = U0;
U0 = U;
dt_1 = dt0;


%% Remaining steps (BDF2)
while t<=Tmax
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

        % remove unused storage
        Tout = Tout(1:jmax);
        Uout = Uout(1:jmax);
        Fout = Fout(1:jmax);
        Xout = Xout(1:jmax);

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


%% Function definitions
% computes G (where G(H,F) = Ht = F - Qx)
function G = calcG(U,F,Dx,Dxx,Dxxx,Dxxxx,R,C,cotb)
U2 = U.*U;
U3 = U2.*U;
U4 = U3.*U;
U5 = U4.*U;
U6 = U5.*U;

DxU = Dx*U;
DxxU = Dxx*U;
DxxxU = Dxxx*U;
DxxxxU = Dxxxx*U;
DxF = Dx*F;

G = F - 2 * U2.*DxU...
    + 2*cotb * U2.*DxU.*DxU...
    + 2*cotb/3 * U3.*DxxU...
    - 1/C * U2.*DxU.*DxxxU...
    - 1/3/C * U3.*DxxxxU...
    - 16*R/5 * U5.*DxU.*DxU...
    - 8*R/15 * U6.*DxxU...
    + 8*R/3 * U3.*DxU.*F...
    + 2*R/3 * U4.*DxF;
end

% computes the Jacobian for G for the Newton solver
function Jg = calcJg(U,F,I,Dx,Dxx,Dxxx,Dxxxx,R,C,cotb)
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

DxU = Dx*U;
DxxU = Dxx*U;
DxxxU = Dxxx*U;
DxxxxU = Dxxxx*U;
DxF = Dx*F;

Jg = -2 * ( J2.*DxU + Dx.*U2 )...
     + 2*cotb * ( J2.*DxU.*DxU + 2*Dx.*DxU.*U2 )...
     + 2*cotb/3 * ( J3.*DxxU + Dxx.*U3 )...
     - 1/C * ( J2.*DxU.*DxxxU + Dx.*DxxxU.*U2 + Dxxx.*DxU.*U2 )...
     - 1/3/C * ( J3.*DxxxxU + Dxxxx.*U3 )...
     - 16*R/5 * ( J5.*DxU.*DxU + 2*Dx.*DxU.*U5 )...
     - 8*R/15 * ( J6.*DxxU + Dxx.*U6 )...
     + 8*R/3 * ( J3.*DxU + Dx.*U3 ) .* F...
     + 2*R/3 * ( J4.*DxF );
end

