function D = derMat(order,N,dx)
%DERMAT constructs derivative matrices
%    returns an NxN matrix for the derivative of
%    the order specified with a uniform grid spacing of dx

v = ones(N,1);

switch order
    case 1
        % central difference, first derivative, [-1/2, 0, 1/2]
        D = spdiags([-v/2,v/2],[-1,1],N,N);
        D(N,1) = 1/2;
        D(1,N) = -1/2;
        D = D/dx;
        return
    case 2
        % central difference, second derivative, [1, -2, 1]
        D = spdiags([v,-2*v,v],[-1,0,1],N,N);
        D(N,1) = 1;
        D(1,N) = 1;
        D = D/(dx^2);
        return
    case 3
        % central difference, third derivative, [-1/2, 1, 0, -1, 1/2]
        D = spdiags([-v/2,v,-v,v/2],[-2,-1,1,2],N,N);
        D(N,1) = -1;
        D(1,N) = 1;
        D(N,2) = 1/2;
        D(2,N) = -1/2;
        D(N-1,1) = 1/2;
        D(1,N-1) = -1/2;
        D = D/(dx^3);
        return
    case 4
        % central difference, fourth derivative, [1, -4, 6, -4, 1]
        D = spdiags([v,-4*v,6*v,-4*v,v],[-2,-1,0,1,2],N,N);
        D(N,1) = -4;
        D(1,N) = -4;
        D(N,2) = 1;
        D(2,N) = 1;
        D(N-1,1) = 1;
        D(1,N-1) = 1;
        D = D/(dx^4);
        return
    case 6
        % central difference, sixth derivative, [1, -6, 15, -20, 15, -6, 1]
        D = spdiags([v,-6*v,15*v,-20*v,15*v,-6*v,v],[-3,-2,-1,0,1,2,3],N,N);
        D(N,1) = 15;
        D(1,N) = 15;
        D(N,2) = -6;
        D(2,N) = -6;
        D(N-1,1) = -6;
        D(1,N-1) = -6;
        D(1,N-2) = 1;
        D(2,N-1) = 1;
        D(3,N) = 1;
        D(N-2,1) = 1;
        D(N-1,2) = 1;
        D(N,3) = 1;
        D = D/(dx^6);
        return
    otherwise
        % order of derivative not supported
        if order < 0
            message = 'order of derivative cannot be negative';
        elseif mod(order,1) > 0
            message = 'order of derivative must be an integer';
        else
            message = ['order of "', num2str(order), '" is not valid'];
        end
        
        error(message);
end


end
