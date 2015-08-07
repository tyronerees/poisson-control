function multdata = genAt(def_setup,multdata,K)
% generates matrices once for geometric multigrid for the advection
% diffusion equation

% Takes the discrete adjoint (the transpose) matrix

for i = 1:def_setup.pow
    multdata(i).Aadj = multdata(i).A';
    multdata(i).GSfwdadj = tril(multdata(i).Aadj,0); % (D-L);
    multdata(i).GSbwdadj = triu(multdata(i).Aadj,0); % (D-U);
end

multdata(def_setup.pow+1).GSfwdadj = tril(K,0); % (D-L);
multdata(def_setup.pow+1).GSbwdadj = triu(K,0); % (D-U);

%     h = 2^(-i+1);           % define h
%     
%     %% get the sd constant, delta
%     hk = min(abs(h/sd.c), abs(h/sd.s));
%     
%     P_h = h*(abs(sd.gamma))/((2)*sd.epsilon);
%     
%     % delta (ESW pg 136)
%     if P_h > 1
%         sd.delta = hk*(1-1/P_h)/((2)*abs(sd.gamma));
%     else
%         sd.delta = 0;
%     end
%     
%     %% assemble the stiffness matrix
%     Abar = stiffsupg(h,sd.c,sd.s,sd.epsilon,sd.gamma,sd.delta);
%     Abar(:,multdata(i).bdn) =[];
%     Abar(multdata(i).bdn,:) = [];
%     Amat(i).A = Abar';
%     
%     %% assemble the smoothing matrices
%     if smoother == 'gs'
% %         D = diag(diag(Amat(i).A)); % define the diagonal
% %         L = tril(-Amat(i).A,-1);   % define the strictly lower part
% %         U = triu(-Amat(i).A,1);    % define the strictly upper part
%         Amat(i).GSfwd = tril(Amat(i).A,0); % (D-L);
%         Amat(i).GSbwd = triu(Amat(i).A,0); % (D-U);
% %         Amat(i).Jfwd = Amat(i).GSfwd\triu(-Amat(i).A,1); % U; save the forward iteration matrix
% %         Amat(i).Jbwd = Amat(i).GSbwd\tril(-Amat(i).A,-1);% L; save the backward iteration matrix
%     end
% end
% 
% if smoother == 'gs'
% %         D = diag(diag(K)); % define the diagonal
% %         L = tril(-K,-1);   % define the strictly lower part
% %         U = triu(-K,1);    % define the strictly upper part
%         Amat(pow+1).GSfwd = tril(K,0); % (D-L);
%         Amat(pow+1).GSbwd = triu(K,0); % (D-U);
% %         Amat(pow+1).Jfwd = Amat(pow+1).GSfwd\triu(-K,1);  % U; save the forward iteration matrix
% %         Amat(pow+1).Jbwd = Amat(pow+1).GSbwd\tril(-K,-1); % L; save the backward iteration matrix
% end
