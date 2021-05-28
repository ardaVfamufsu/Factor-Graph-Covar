function [delta_x,rhs_x,rhs_z,w_A,unw_A,W,unw_R] =...
    linearSLAM_newWeighting(x_vec,Y,u,n,G,F,x0,Se_sqr,Sp_sqr)

%% Inputs:
%%  x_vec0 is the initial state (linearization point) should be size 4xn
%%  Y are the measurements (size 2xn)
%%  control_mat (size 2xn) inputs (u)
%%  G -- Dynamics matrix
%%  F -- 2x4 measurement matrix
%%  x0 -- state at time 0 (before time 1), assumed fixed/known
%%  Se_sqr -- Sqrt of Measurement noise covariance
%%  Sp_sqr -- Sqrt of dynamics noise covariance
%% Returns:
%%  del_x_sol -- delta x.  How much to move things by (so is it really linear?)
%%  rhs_x -- the right hand side, weighted part, corresponding to dynamics
%%  rhs_z -- the right hadn size, weighted, part corresponding to measurements
%%  w_A -- The _weighted_ A matrix
%%  unw_A -- the A matrix, _unweighted_
%%  W -- The weight matrix.
%%  unw_R -- the unweighted residual matrix
%%
%% In general, delta_x = inv(A^T W A) A^T W unw_R


x1=reshape(x_vec(1:4*n),[4,n]);
%% First, create the "right hand size" (the weighted and unweighted residual vectors)
unw_rhs_x = zeros(n*4,1);
rhs_x = zeros(n*4,1);
% rhs vector for x
for i=1:n
    if i==1
        a_vec_i=x1(:,i)-G*x0-u(:,i);
    else
        a_vec_i=x1(:,i)-G*x1(:,i-1)-u(:,i);
    end
    unw_rhs_x((i-1)*4+1:i*4) = a_vec_i;
    rhs_x(i*4-3:i*4) = Sp_sqr*a_vec_i;
end

% rhs vector for z
unw_rhs_z = zeros(n*2,1);
rhs_z = zeros(n*2,1);
for i=1:n
  v_states=x1(:,i);
  v_meas=Y(:,i);
  c_vec_i=v_meas-F*v_states;  
  unw_rhs_z(i*2-1:i*2)=c_vec_i;
  rhs_z(i*2-1:i*2) = Se_sqr*c_vec_i;
end
% combine the rhs vectors
rhs=[rhs_x; rhs_z];
unw_R = [unw_rhs_x; unw_rhs_z];

%% Second, create the A and W matrices
% A matrix for x
Ax_mid=[G -eye(4)];
for i=1:n
    if i==1
        Ax=[-eye(4) repmat(zeros(4), 1,n-1)];
    else
        Ax=[Ax; repmat(zeros(4), 1,i-2)  Ax_mid repmat(zeros(4), 1,n-i)];
    end    
end
% A matrix for z
Az=[];
for i=1:n
  Az1=[repmat(zeros(2,4), 1,i-1)  F repmat(zeros(2,4), 1,n-i)];
  Az=[Az; Az1];  
end

% combine the A matrices
unw_A=[Ax; Az];

% Generate the W matrix
% Making sparse as it is mostly zeros
n_Sp = length(Sp_sqr(:));
n_Se = length(Se_sqr(:));
r_Sp = size(Sp_sqr,1);
c_Sp = size(Sp_sqr,2);
r_Se = size(Se_sqr,1);
c_Se = size(Se_sqr,2);
nnz_entries = n_Sp*n + n_Se*n;
data_arr = zeros(nnz_entries,1);
row_arr = zeros(nnz_entries,1);
col_arr = zeros(nnz_entries,1);
% Do the Sp entries
for ii=1:n
    [rows,cols] = create_ij(ii*4-3,ii*4-3,r_Sp,c_Sp);
    data_arr(ii*n_Sp-n_Sp+1:ii*n_Sp) = Sp_sqr(:);
    row_arr(ii*n_Sp-n_Sp+1:ii*n_Sp) = rows(:);
    col_arr(ii*n_Sp-n_Sp+1:ii*n_Sp) = cols(:);
end
% And the Se entries
zmi = n*4; %Z's matrix index
zai = n_Sp*n; %Z's array index
for ii=1:n
    [rows,cols] = create_ij(zmi+ ii*2-1, zmi+ ii*2-1, r_Se,c_Se);
    data_arr(zai + ii*n_Se-n_Se+1: zai+ii*n_Se) = Se_sqr(:);
    row_arr(zai + ii*n_Se-n_Se+1: zai+ii*n_Se) = rows(:);
    col_arr(zai + ii*n_Se-n_Se+1: zai+ii*n_Se) = cols(:);
end
W_sqr = sparse(row_arr,col_arr,data_arr);
W = W_sqr.' * W_sqr;
   
%% Solve for delta_x
w_A = W_sqr*unw_A;
AtA = w_A.'*w_A;
Aty = w_A.'*rhs;
delta_x = AtA\Aty;

function [rows,cols] = create_ij(row_offset,col_offset, n_rows,n_cols)
%% Create an indexing for a matrix of size n_rows,n_cols where the top
%% left corner is at row_offset,col_offset
    [local_c,local_r] = meshgrid(0:n_cols-1, 0:n_rows-1);
    rows = row_offset+local_r(:);
    cols = col_offset+local_c(:);

