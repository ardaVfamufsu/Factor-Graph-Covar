function [x_sol, rhs_x, rhs_z, Ax, Az, A,rhs, Aw,rhs_w,W,rhs_x_unw,rhs_z_unw,R,Apw]=...
    linearSLAM_noLandmark_gps_Mestimate(x_vec,Y,u,n,G,F,x0,Se_sq,Sp_sq)
% x_vec: linearization point for the state vector
% x0: vehicle states at time 0

%vehicle trajectory, 4 by n matrix
x1=reshape(x_vec(1:4*n),[4,n]); 
% rhs vector for x
for i=1:n
    if i==1
        a_vec_i=x1(:,i)-G*x0-u(:,i);
        rhs_x=Sp_sq*a_vec_i; % weighted right hand side
        rhs_x_unw=a_vec_i;   % unweighted right hand side, needed for M estimator (prob not needed)
    else
        a_vec_i=x1(:,i)-G*x1(:,i-1)-u(:,i);
        rhs_x=[rhs_x; Sp_sq*a_vec_i];
        rhs_x_unw=[rhs_x_unw; a_vec_i];
    end
end
% rhs vector for z
rhs_z=[];
% unweighted right hand side, needed for M estimator
rhs_z_unw=[];
for i=1:n
  v_states=x1(:,i);
  v_meas=Y(:,i);
  c_vec_i=v_meas-F*v_states;  
  rhs_z=[rhs_z; Se_sq*c_vec_i];  
  rhs_z_unw=[rhs_z_unw; c_vec_i];  
end
% combine the rhs vectors
rhs=[rhs_x; rhs_z];
rhs_unw=[rhs_x_unw; rhs_z_unw];
% A matrix for x
Ax_mid=[Sp_sq*G -Sp_sq*eye(4)];
for i=1:n
    if i==1
        Ax=[-Sp_sq*eye(4) repmat(zeros(4), 1,n-1)];
    else
        Ax=[Ax; repmat(zeros(4), 1,i-2)  Ax_mid repmat(zeros(4), 1,n-i)];
    end    
end
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

% robust estimate of scale- (Montgomery et al., 2006, p. 374)
% apply to residuals of measurements only
s_robust=median(abs(rhs_z-median(rhs_z)))/0.6745;
z_res=rhs_z/s_robust;

% Ramsay's Ea function with a=0.3
%Ra=0.3;
%wHub_z=exp(-Ra*abs(z_res));
% % Cauchy function
k_c=1.645;
wHub_z=1./(1+(z_res/k_c).^2);
% % Tukey bisquare
% k_bsq=3;
% for i=1:length(z_res)
%     if abs(z_res(i))<=k_bsq
%         wHub_z(i,1)=(1-z_res(i)^2/k_bsq^2)^2;
%     else
%         wHub_z(i,1)=1e-6;
%     end
% end
weights=[ones(4*n,1);wHub_z];
%weights=[ones(4*n,1);ones(length(z_res),1)]; % unweighted
% find the weights matrix
W=diag(weights);
Wsq=diag(sqrt(weights));
% find weighted A matrix and rhs vector 
Aw=Wsq*A;
rhs_w=Wsq*rhs;

% reorder columns of A to exploit sparsity
indcolA=colamd(Aw);
Apw=Aw(:,indcolA);
% inverse reordering to find original variables
[~,indInvSortA]=sort(indcolA);
        
% solve for the states
if rank(Apw)==size(Apw,2)
    % upper triangular matrix R obtained by cholesky factorization of Ap
    R=chol(Apw'*Apw);
    % solve for x using back substitution
    xp = R\(R'\(Apw'*rhs_w));
else
    %QR factorization 
    xp=Apw\rhs_w; 
end

% % use LDL decompositon, alternatively
% %[L,D]=ldl(Ap'*Ap);
% % solve for x using back substitution
% %xp = (D*L')\(L\(Ap'*rhs));
x_sol_factor=xp(indInvSortA);
x_sol=x_sol_factor;


end


