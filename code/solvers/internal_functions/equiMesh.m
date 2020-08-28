function [Z,P,dxmin] = equiMesh(X,W,L,K,d,dxinf,dxsup)
%EQUIMESH Creates a locally uniform mesh equidistributing a weighting W
%   EQUIMESH uses the method described in the relevant PDF to produce a
%            sub-equidistributing mesh for W.
%
%   INPUTS:  X - point spacing for W. Should be a strictly increasing sequence.
%            W - weighting function (to be equidistributed). Must be
%                non-negative.
%            L - length of domain. Should be greater than X(end)
%            K - maximum ratio between adjacent spaces
%            d - value of W at which dxmin will be 90% of the way to dxinf.
%                If d = 0 then dxmin = dxmax
%        dxinf - limit of dx as W -> infinity
%        dxsup - value of dx at W = 0
%
%
%   OUTPUTS: Z - output mesh
%            P - padded version of W (on the same grid W)
%        dxmin - minimum value of dx

% check weighting function is allowed
if sum(W<0)>0
    error('Weighting function must be non-negative')
end

Nmin = floor(L/dxinf);
if sum(W)==0
    Z = linspace(0,L,Nmin+1)';
    Z = Z(1:end-1);
    P = zeros(size(Z));
    dxmin = Z(2)-Z(1);
    return
end

% if dxinf = dxsup the the grid should be uniform
if dxsup==dxinf
    Z = linspace(0,L,Nmin+1)';
    Z = Z(1:end-1);
    P = zeros(size(Z));
    dxmin = Z(2)-Z(1);
    return
end


%% Prep inputs
N = length(X);
wmax = max(W);

if d==0
    dxmin = dxinf;
else
    d = tan(0.9)/d;
    dxmin = dxsup + (dxinf-dxsup)*tanh(d*wmax);
end


a = (dxmin*wmax)/(dxsup-dxmin); % additive damping paramter
c = wmax*dxmin*dxsup/(dxsup-dxmin); %

% account for periodicity
X3 = [X-L;X;X+L];
W3 = [W;W;W] + a;

% convert to a MATLAB-style mesh
O = ones(3*N,3*N);
W3mesh = O.*W3;
X3mesh = O.*X3;


%% Create grid
l = log(K)/c; % padding parameter
G3 = W3mesh'./(1+l*abs(X3mesh'-X3mesh).*W3mesh'); % family of envelope functions
P3 = max(G3,[],2); % padded weight function (non-periodic)
P = max([[P3(1:N);P3(N+1)],[P3(N+1:2*N);P3(2*N+1)],[P3(2*N+1:3*N);P3(1)]],[],2); % wrap to periodic
Q = cumtrapz([X;L],P);

n = ceil(Q(end)/c); % number of gridpoints

%n = max(n,64);
%n = min(n,512);

Y = linspace(0,Q(end),n+1);

Z = interp1(Q,[X;L],Y,'pchip')';
Z = Z(1:end-1);

P = P(1:end-1);
end

