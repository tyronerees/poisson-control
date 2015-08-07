function [dir,neu] = robbdy(N)

k = 1:N;
bdy1 = [k' k'+1];
bdy2 = [k'*(N+1) (k'+1)*(N+1)];
bdy3 = [(N+1)^2-k'+1 (N+1)^2-k'];
k = N:-1:1;
bdy4 = [k'*(N+1)+1 (k'-1)*(N+1)+1];                 % define bdyn, which gives you the node numbers
                                                    % touching the nth edge of the boundary

dir = [bdy4;bdy1];                                  % Define Dirichlet bc's on edge 1&2
neu = [bdy2;bdy3];                                  % Define Neumann bc's on edge 3&4