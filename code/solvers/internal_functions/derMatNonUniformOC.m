function D = derMatNonUniformOC(ord,acc,X,L)
%DERMATNONUNIFORMOC constructs non-centred derivative matrices
%    returns an NxN matrix for the derivative of
%    the order and accuracy specified on the non-uniform grid X.

H0 = X - circshift(X,1);
H0(1) = H0(1) + L;

H_1 = circshift(H0,1);
H1 = circshift(H0,-1);
H2 = circshift(H0,-2);
H3 = circshift(H0,-3);

N = length(X);

switch ord
    case 1
        switch acc
            case 2
                v_1 = -H1./H0./(H0+H1);
                v0 = (H1-H0)./H0./H1;
                v1 = H0./H1./(H0+H1);

                % probably can remove these shifts by cancelling with an earlier one
                v1 = circshift(v1,1);
                v_1 = circshift(v_1,-1);

                % diagonals
                D = spdiags([v_1,v0,v1],[-1,0,1],N,N);

                % corners
                D(1,N) = v_1(N);
                D(N,1) = v1(1);
                return
            otherwise
                message = ['accuracy of "', num2str(acc), '" is not valid for order ', num2str(ord)];
                error(message);
        end
    case 2
        switch acc
            case 1
%                 v_1 = 2./H0./(H0 + H1);
%                 v0 = -2./H0./H1;
%                 v1 = 2./H1./(H0 + H1);
%                 
%                 % probably can remove these shifts by cancelling with an earlier one
%                 v1 = circshift(v1,1);
%                 v_1 = circshift(v_1,-1);
%                 
%                 % diagonals
%                 D = spdiags([v_1,v0,v1],[-1,0,1],N,N);
%                 
%                 % corners
%                 D(1,N) = v_1(N);
%                 D(N,1) = v1(1);
%                 return
            case 2
                v_1 = (2*(2*H1 + H2))./(H0.*(H0 + H1).*(H0 + H1 + H2));
                v0 = -(2*(2*H1 - H0 + H2))./(H0.*H1.*(H1 + H2));
                v1 = (2*(H1 - H0 + H2))./(H1.*H2.*(H0 + H1));
                v2 = (2*(H0 - H1))./(H2.*(H1 + H2).*(H0 + H1 + H2));

                % probably can remove these shifts by cancelling with an earlier one
                v2 = circshift(v2,2);
                v1 = circshift(v1,1);
                v_1 = circshift(v_1,-1);

                % diagonals
                D = spdiags([v_1,v0,v1,v2],[-1,0,1,2],N,N);

                % corners
                D(1,N) = v_1(N);

                D(N,1) = v1(1);
                D(N,2) = v2(2);
                D(N-1,1) = v2(1);
                return
%             case 3
%                 v_2 = (2*(H1.^2 - 2*H0.*H1 - H0.*H2 + H1.*H2))./(H_1.*(H0 + H_1).*(H0 + H1 + H_1).*(H0 + H1 + H2 + H_1));
%                 v_1 = (2*(2*H0.*H1 - H1.^2 + H0.*H2 - H1.*H2 + 2*H1.*H_1 + H2.*H_1))./(H0.*H_1.*(H0 + H1).*(H0 + H1 + H2));
%                 v0 = (2*(H0.^2 + H1.^2 - 4*H0.*H1 - 2*H0.*H2 + H1.*H2 + H0.*H_1 - 2*H1.*H_1 - H2.*H_1))./(H0.*H1.*(H1 + H2).*(H0 + H_1));
%                 v1 = (2*(2*H0.*H1 - H0.^2 + 2*H0.*H2 - H0.*H_1 + H1.*H_1 + H2.*H_1))./(H1.*H2.*(H0 + H1).*(H0 + H1 + H_1));
%                 v2 = (2*(H0.^2 - 2*H0.*H1 + H0.*H_1 - H1.*H_1))./(H2.*(H1 + H2).*(H0 + H1 + H2).*(H0 + H1 + H2 + H_1));
% 
%                 % probably can remove these shifts by cancelling with an earlier one
%                 v2 = circshift(v2,2);
%                 v1 = circshift(v1,1);
%                 v_1 = circshift(v_1,-1);
%                 v_2 = circshift(v_2,-2);
% 
%                 % diagonals
%                 D = spdiags([v_2,v_1,v0,v1,v2],[-2,-1,0,1,2],N,N);
% 
%                 % corners
%                 D(1,N) = v_1(N);
%                 D(1,N-1) = v_2(N-1);
%                 D(2,N) = v_2(N);
% 
%                 D(N,1) = v1(1);
%                 D(N,2) = v2(2);
%                 D(N-1,1) = v2(1);
%                 return
            otherwise
                message = ['accuracy of "', num2str(acc), '" is not valid for order ', num2str(ord)];
                error(message);
        end
    case 3
        switch acc
            case 2
                v_1 = -(6*(3*H1 + 2*H2 + H3))./(H0.*(H0 + H1).*(H0 + H1 + H2).*(H0 + H1 + H2 + H3));
                v0 = (6*(3*H1 - H0 + 2*H2 + H3))./(H0.*H1.*(H1 + H2).*(H1 + H2 + H3));
                v1 = -(6*(2*H1 - H0 + 2*H2 + H3))./(H1.*H2.*(H0 + H1).*(H2 + H3));
                v2 = (6*(2*H1 - H0 + H2 + H3))./(H2.*H3.*(H1 + H2).*(H0 + H1 + H2));
                v3 = -(6*(2*H1 - H0 + H2))./(H3.*(H2 + H3).*(H1 + H2 + H3).*(H0 + H1 + H2 + H3));

                % probably can remove these shifts by cancelling with an earlier one
                v3 = circshift(v3,3);
                v2 = circshift(v2,2);
                v1 = circshift(v1,1);
                v_1 = circshift(v_1,-1);

                % diagonals
                D = spdiags([v_1,v0,v1,v2,v3],[-1,0,1,2,3],N,N);

                % corners
                D(1,N) = v_1(N);

                D(N,1) = v1(1);
                D(N,2) = v2(2);
                D(N-1,1) = v2(1);
                D(N,3) = v3(3);
                D(N-1,2) = v3(2);
                D(N-2,1) = v3(1);
                return
            otherwise
                message = ['accuracy of "', num2str(acc), '" is not valid for order ', num2str(ord)];
                error(message);
        end
    case 4
        switch acc
%             case 1
%                 v_2 = 24./(H_1.*(H0 + H_1).*(H0 + H1 + H_1).*(H0 + H1 + H2 + H_1));
%                 v_1 = -24./(H0.*H_1.*(H0 + H1).*(H0 + H1 + H2));
%                 v0 = 24./(H0.*H1.*(H1 + H2).*(H0 + H_1));
%                 v1 = -24./(H1.*H2.*(H0 + H1).*(H0 + H1 + H_1));
%                 v2 = 24./(H2.*(H1 + H2).*(H0 + H1 + H2).*(H0 + H1 + H2 + H_1));
% 
%                 % probably can remove these shifts by cancelling with an earlier one
%                 v2 = circshift(v2,2);
%                 v1 = circshift(v1,1);
%                 v_1 = circshift(v_1,-1);
%                 v_2 = circshift(v_2,-2);
% 
%                 % diagonals
%                 D = spdiags([v_2,v_1,v0,v1,v2],[-2,-1,0,1,2],N,N);
% 
%                 % corners
%                 D(1,N) = v_1(N);
%                 D(1,N-1) = v_2(N-1);
%                 D(2,N) = v_2(N);
% 
%                 D(N,1) = v1(1);
%                 D(N,2) = v2(2);
%                 D(N-1,1) = v2(1);
%                 return
%             case 2
%                 v_3 = (24*(3*H1 - 2*H0 + 2*H2 + H3 - H_1))./(H_2.*(H_1 + H_2).*(H0 + H1 + H_1 + H_2).*(H0 + H1 + H2 + H_1 + H_2).*(H0 + H1 + H2 + H3 + H_1 + H_2));
%                 v_2 = (24*(2*H0 - 3*H1 - 2*H2 - H3 + H_1 + H_2))./(H_1.*H_2.*(H0 + H1 + H_1).*(H0 + H1 + H2 + H_1).*(H0 + H1 + H2 + H3 + H_1));
%                 v_1 = -(24*(2*H0 - 3*H1 - 2*H2 - H3 + 2*H_1 + H_2))./(H_1.*(H0 + H1).*(H_1 + H_2).*(H0 + H1 + H2).*(H0 + H1 + H2 + H3));
%                 v1 = (24*(3*H0 - 2*H1 - 2*H2 - H3 + 2*H_1 + H_2))./(H2.*(H0 + H1).*(H2 + H3).*(H0 + H1 + H_1).*(H0 + H1 + H_1 + H_2));
%                 v2 = -(24*(3*H0 - 2*H1 - H2 - H3 + 2*H_1 + H_2))./(H2.*H3.*(H0 + H1 + H2).*(H0 + H1 + H2 + H_1).*(H0 + H1 + H2 + H_1 + H_2));
%                 v3 = (24*(3*H0 - 2*H1 - H2 + 2*H_1 + H_2))./(H3.*(H2 + H3).*(H0 + H1 + H2 + H3).*(H0 + H1 + H2 + H3 + H_1).*(H0 + H1 + H2 + H3 + H_1 + H_2));
% 
%                 % probably can remove these shifts by cancelling with an earlier one
%                 v3 = circshift(v3,3);
%                 v2 = circshift(v2,2);
%                 v1 = circshift(v1,1);
%                 v_1 = circshift(v_1,-1);
%                 v_2 = circshift(v_2,-2);
%                 v_3 = circshift(v_3,-3);
% 
%                 % diagonals
%                 D = spdiags([v_3,v_2,v_1,v1,v2,v3],[-3,-2,-1,1,2,3],N,N);
% 
%                 % corners
%                 D(1,N) = v_1(N);
%                 D(1,N-1) = v_2(N-1);
%                 D(2,N) = v_2(N);
%                 D(1,N-2) = v_3(N-2);
%                 D(2,N-1) = v_3(N-1);
%                 D(3,N) = v_3(N);
% 
%                 D(N,1) = v1(1);
%                 D(N,2) = v2(2);
%                 D(N-1,1) = v2(1);
%                 D(N,3) = v3(3);
%                 D(N-1,2) = v3(2);
%                 D(N-2,1) = v3(1);
%                 return
            case 2
                v_2 = (24*(3*H1 - H0 + 2*H2 + H3))./(H_1.*(H0 + H_1).*(H0 + H1 + H_1).*(H0 + H1 + H2 + H_1).*(H0 + H1 + H2 + H3 + H_1));
                v_1 = -(24*(3*H1 - H0 + 2*H2 + H3 - H_1))./(H0.*H_1.*(H0 + H1).*(H0 + H1 + H2).*(H0 + H1 + H2 + H3));
                v0 = (24*(3*H1 - 2*H0 + 2*H2 + H3 - H_1))./(H0.*H1.*(H1 + H2).*(H0 + H_1).*(H1 + H2 + H3));
                v1 = -(24*(2*H1 - 2*H0 + 2*H2 + H3 - H_1))./(H1.*H2.*(H0 + H1).*(H2 + H3).*(H0 + H1 + H_1));
                v2 = (24*(2*H1 - 2*H0 + H2 + H3 - H_1))./(H2.*H3.*(H1 + H2).*(H0 + H1 + H2).*(H0 + H1 + H2 + H_1));
                v3 = (24*(2*H0 - 2*H1 - H2 + H_1))./(H3.*(H2 + H3).*(H1 + H2 + H3).*(H0 + H1 + H2 + H3).*(H0 + H1 + H2 + H3 + H_1));
                
                % probably can remove these shifts by cancelling with an earlier one
                v3 = circshift(v3,3);
                v2 = circshift(v2,2);
                v1 = circshift(v1,1);
                v_1 = circshift(v_1,-1);
                v_2 = circshift(v_2,-2);

                % diagonals
                D = spdiags([v_2,v_1,v0,v1,v2,v3],[-2,-1,0,1,2,3],N,N);

                % corners
                D(1,N) = v_1(N);
                D(1,N-1) = v_2(N-1);
                D(2,N) = v_2(N);

                D(N,1) = v1(1);
                D(N,2) = v2(2);
                D(N-1,1) = v2(1);
                D(N,3) = v3(3);
                D(N-1,2) = v3(2);
                D(N-2,1) = v3(1);
                return
            otherwise
                message = ['accuracy of "', num2str(acc), '" is not valid for order ', num2str(ord)];
                error(message);
        end
    otherwise
        % order of derivative not supported
        if ord < 0
            message = 'order of derivative cannot be negative';
        elseif mod(ord,1) > 0
            message = 'order of derivative must be an integer';
        else
            message = ['order of "', num2str(ord), '" is not supoorted'];
        end
        error(message);
end

end


