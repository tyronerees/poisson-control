function [dir,neu] = neubdy(N)

k = 1:N;
bdy1 = [k' k'+1];
bdy2 = [k'*(N+1) (k'+1)*(N+1)];
bdy3 = [(N+1)^2-k'+1 (N+1)^2-k'];
k = N:-1:1;
bdy4 = [k'*(N+1)+1 (k'-1)*(N+1)+1];                 % define bdyn, which gives you the node numbers
                                                    % touching the nth edge of the boundary

%  keyboard
neu = [bdy1;bdy2(1:end-1,:);bdy3(2:end,:);bdy4];             % Define Neumann bc's everywhere except on one element
dir = bdy3(1,1);                                  % Make this element Dirichlet                