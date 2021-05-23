function [x_sol, rhs_x, rhs_z, Ax, Az, A,rhs,R,Ap]=...
  linearSLAM_forwardVel_noLandmark_gps(x_vec,Y,u,n,F,x0,Se_sq,Sp_sq,T)
% x_vec: linearization point for the state vector
% x0: vehicle states at time 0

%vehicle trajectory, 4 by n matrix
x1=reshape(x_vec(1:4*n),[4,n]); 
% rhs vector for x
for i=1:n
    if i==1
        G=[ 1 0 -x0(4)*T*sin(x0(3)) T*cos(x0(3));
            0 1 x0(4)*T*cos(x0(3)) T*sin(x0(3));
            0 0 1 0;
            0 0 0 1];
        %a_vec_i=x1(:,i)-G*x0-u(:,i);
        f_prev=[x0(1)+x0(4)*T*cos(x0(3));
                x0(2)+x0(4)*T*sin(x0(3));
                x0(3);
                x0(4)];
        a_vec_i=x1(:,i)-f_prev-u(:,i);
        rhs_x=Sp_sq*a_vec_i;
        Ax=[-Sp_sq*eye(4) repmat(zeros(4), 1,n-1)];
    else
        G=[ 1 0 -x1(4,i-1)*T*sin(x1(3,i-1)) T*cos(x1(3,i-1));
        0 1 x1(4)*T*cos(x1(3,i-1)) T*sin(x1(3,i-1));
        0 0 1 0;
        0 0 0 1];
        %a_vec_i=x1(:,i)-G*x1(:,i-1)-u(:,i);
        f_prev=[x1(1,i-1)+x1(4,i-1)*T*cos(x1(3,i-1));
                x1(2,i-1)+x1(4,i-1)*T*sin(x1(3,i-1));
                x1(3,i-1);
                x1(4,i-1)];
        a_vec_i=x1(:,i)-f_prev-u(:,i);
        rhs_x=[rhs_x; Sp_sq*a_vec_i];
        Ax_mid=[Sp_sq*G -Sp_sq*eye(4)];
        Ax=[Ax; repmat(zeros(4), 1,i-2)  Ax_mid repmat(zeros(4), 1,n-i)]; 
    end
end
% rhs vector for z
rhs_z=[]; 
for i=1:n
  v_states=x1(:,i);
  v_meas=Y(:,i);
  c_vec_i=v_meas-F*v_states;  
  rhs_z=[rhs_z; Se_sq*c_vec_i];  
end
% combine the rhs vectors
rhs=[rhs_x; rhs_z];

% A matrix for x
%Ax_mid=[Sp_sq*G -Sp_sq*eye(4)];
%for i=1:n
%    if i==1
%        Ax=[-Sp_sq*eye(4) repmat(zeros(4), 1,n-1)];
%    else
%        Ax=[Ax; repmat(zeros(4), 1,i-2)  Ax_mid repmat(zeros(4), 1,n-i)];
%    end    
%end
% A matrix for z
Az=[];
for i=1:n
  Az1=[repmat(zeros(2,4), 1,i-1)  Se_sq*F repmat(zeros(2,4), 1,n-i)];
  Az=[Az; Az1];  
end

% combine the A matrices
A=[Ax; Az];
% solve for the states
%x_sol=A\rhs;

%% factor graph approach
% reorder columns of A to exploit sparsity
indcolA=colamd(A);
Ap=A(:,indcolA);
% inverse reordering to find original variables
[~,indInvSortA]=sort(indcolA);

if rank(A)==size(A,2)
    % upper triangular matrix R obtained by cholesky factorization of Ap
    R=chol(Ap'*Ap);
    % solve for x using back substitution
    xp = R\(R'\(Ap'*rhs));
else
    %QR factorization 
    %xp=Ap\rhs; 
    [~,R] = qr(Ap);
    xp = R\(R'\(Ap'*rhs));
end

% % use LDL decompositon, alternatively
% %[L,D]=ldl(Ap'*Ap);
% % solve for x using back substitution
% %xp = (D*L')\(L\(Ap'*rhs));
x_sol_factor=xp(indInvSortA);
x_sol=x_sol_factor;
% % reorder columns and rows of R matrix
% R_sol1=R(:,indInvSortA);
% R_sol=R_sol1(indInvSortA,:);


end


