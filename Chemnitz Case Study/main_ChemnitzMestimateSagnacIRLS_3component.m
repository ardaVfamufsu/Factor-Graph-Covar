function [x_sol, rhs_x, rhs_z, Ax, Az, A,rhs,sig_pHat1,sig_pHat2, sig_eHat]=main_ChemnitzMestimateSagnacIRLS_3component(flag_unbiased,flag_scale_cov_change)
% fl_unbiased: 1 to use unbiased variance estimator, 0 to use biased variance estimator
% flag_scale_cov_change:1 if trying to "scale" adjusted covariance to keep with order of magnitude change
% input data
%M=readmatrix('Data_Chemnitz'); % columns are: time stamp, pseudorange, x,y,z coordinates of satellite
%prange_mat=M(1:90273,[2 3 5 6 7]); 
%gtruth_mat=M(98844:end,[2 3 4 5]);
M=csvread('Data_Chemnitz.csv',0,1);
prange_mat=M(1:90273,[1 2 4 5 6]);% columns are: time stamp, pseudorange, x,y,z coordinates of satellite
gtruth_mat=M(98844:end,[1 2 3 4]);% columns are: time stamp, x,y,z coordinates of vehicle
% get portions of interest
tmeas=prange_mat(:,1);
rmeas=prange_mat(:,2);
Scoord=prange_mat(:,3:5);
coord_true=gtruth_mat(:,2:4);
% measurement instances, no of satellites
ttab=tabulate(tmeas);
t_unique=ttab(:,1);
n_sat=ttab(:,2);
n=length(t_unique);
% as a test
% n=2000;
% state vector 5 by 1: x,y,z of vehicle, b(clock error), d (clock drift)
n_state=5;
% scaling factor for time
scf=1e9;% convert c to m/ns
% speed of light (m/s)
c=3e8;
c=c/scf; 
% receiver clock state - process error covariance (Bar Shalom, p . 511)
Sb=2*2e-19;
Sd=8*pi^2*2e-20;
Sb=Sb*scf^2;% convert to ms or ns
Sd=Sd*scf^2;% convert to ms or ns
% pseudoranges, satellite coordinates
j_all=0;
for i=1:n   
    j=0;
    while j<n_sat(i)
        j=j+1;
        j_all=j_all+1;
        Yn(j,i)=rmeas(j_all);  
        x_sat(j,i)=Scoord(j_all,1);
        y_sat(j,i)=Scoord(j_all,2);
        z_sat(j,i)=Scoord(j_all,3);
        tm(j,i)=tmeas(j_all);
    end
end 
% clear large variables that will not be further used
clear M prange_mat Scoord
% known initial state (one time period before the measurements start)
%   x y z b d  (assume initial clock error b=0. x, y, z are immaterial, because process model does not contain them.)
x0=[0 0 0 0 0]'; 
% iteratively reweighted LS (Montgomery et al., 2006, p. 374)
% iteration parameters
err_sigSol=99;
c1=1;
% linearization point
x_vec0=repmat(zeros(n_state,1),n,1);
sig_eHat=10^2;
% for scf=1e9
%sig_pHat1=Sb*1e-6; %offset error (gives negative variance)
sig_pHat1=Sb; %offset error 
sig_pHat2=Sd*0.01; %drift error

% % for scf=1 (gives negative variances)
% sig_pHat1=1e-10; %offset error
% sig_pHat2=1e-10; %drift error

% remove zero rows from Ax and rhs_x, avoids rank deficiency
cind=1;
for i=1:n_state*(n-1)  % end at n-1 since we omit first n_state data from rhs_x
    if ismember(mod(i,5),[4 0])
        ind_x_keep(cind)=i;
        cind=cind+1;
    end
end
%ind_x_keep=1:n_state*(n-1); % do not remove, causes rank deficiency
mx=length(ind_x_keep);
m_iter=15;
m2_iter=10;
conv_criteria=5;
while c1<m2_iter &&err_sigSol>conv_criteria
    Qc=[sig_pHat1 0; 0 sig_pHat2];
    Q_sq_hat=(inv(Qc))^0.5;
    Sp_sq_hat=[zeros(3,5);
                zeros(2,3) Q_sq_hat];      
    x_vec0_prev=x_vec0;        
    for i=1:m_iter        
        [del_x_sol, rhs_x, rhs_z, Ax, Az, A,rhs,Apw,rhs_w,w_vec,~]=getGPSgraph(x_vec0,Yn,n,x0,sig_eHat,...
            Sp_sq_hat,n_sat,n_state,x_sat,y_sat,z_sat,c,tm);
        x_vec0=x_vec0+del_x_sol;
        %i
        disp(['state estimation, with fixed variance=',num2str(dot(del_x_sol,del_x_sol))])    
    end
    % indices for rhs_z
    m=length(rhs);
    ind_z=n_state*(n-1)+1:m;
    mz=length(ind_z);
    % residuals
    rhs_x1_w=rhs_w(ind_x_keep(1:2:end));
    rhs_x2_w=rhs_w(ind_x_keep(2:2:end));
    rhs_z_w=rhs_w(ind_z);
    % update previous values
    sig_pHat1_prev=sig_pHat1;
    sig_pHat2_prev=sig_pHat2;
    sig_eHat_prev=sig_eHat;
    if flag_unbiased
        disp('unbiased variance estimation - may take ~10 minutes')
        % find new standardized variances - unbiased 
        Apw_x1=Apw(ind_x_keep(1:2:end),:);
        Apw_x2=Apw(ind_x_keep(2:2:end),:);
        Apw_z=Apw(ind_z,:);       
        % using efficient approach
        [sig_pHat1_std, sig_pHat2_std, sig_eHat_std]=getUnbVarEstSLAM(rhs_x1_w,rhs_x2_w,rhs_z_w,Apw_x1,Apw_x2,Apw_z,flag_scale_cov_change);               
    else
        disp('variance estimation')
        % biased (robust) estimate of standardized measurement, process and combined variances
        mx1=length(rhs_x1_w); mx2=length(rhs_x2_w);
        sig_pHat1_std=1/(mx1-1)*dot(rhs_x1_w,rhs_x1_w);  
        sig_pHat2_std=1/(mx2-1)*dot(rhs_x2_w,rhs_x2_w);  
        sig_eHat_std=1/(mz-1)*dot(rhs_z_w,rhs_z_w);
    end
    sig_eHat=(sig_eHat_std)*sig_eHat_prev;
    sig_pHat1=(sig_pHat1_std)*sig_pHat1_prev;
    sig_pHat2=(sig_pHat2_std)*sig_pHat2_prev; 

    % change in variance estimates
    err_vec=[sig_eHat-sig_eHat_prev; sig_pHat1-sig_pHat1_prev; sig_pHat2-sig_pHat2_prev];
    err_sigSol=dot(err_vec,err_vec);  
    % change in state estimates
    err_xSol=dot(x_vec0-x_vec0_prev,x_vec0-x_vec0_prev);
    disp(['end of state & variance estimation, iteration ',num2str(c1)]);
    disp('process and measure variances')
    disp([sig_pHat1, sig_pHat2, sig_eHat])
    disp('process and measure standardized variances')
    disp([sig_pHat1_std,sig_pHat2_std, sig_eHat_std])    
    disp('variance convergence criterion')    
    disp(err_sigSol)
    disp('state convergence criterion')    
    disp(err_xSol)
    % iteration counter
    c1=c1+1;    
end

%vehicle trajectory
x_sol=x_vec0;
v_est=reshape(x_sol,[n_state,n]);

figure(1);
x_origin=coord_true(1,1);
y_origin=coord_true(1,2);
plot(coord_true(:,1)-x_origin,coord_true(:,2)-y_origin);
hold on;
plot(v_est(1,:)-x_origin,v_est(2,:)-y_origin,'.r');
xlabel('x coord')
ylabel('y coord')
legend('Ground Truth','Estimated')
grid on
pbaspect([1 1 1])
set(gca,'LooseInset',get(gca,'TightInset'));
%saveas(gcf,'psrangeClSagnacMestIRLS.jpg');

disp('position estimation error (m), mean, median, 97.5% and max')
e_x=coord_true(1:n,1)-v_est(1,:)';
e_y=coord_true(1:n,2)-v_est(2,:)';
e_dist=sqrt(e_x.^2+e_y.^2);
disp(mean(e_dist));
disp(prctile(e_dist,[50, 97.5]));
disp(max(e_dist))

% errors, residuals, weights
figure(3);
plot(abs(rhs_z_w), w_vec,'o');
xlabel('absolute weighted res')
ylabel('weight')
figure(4);
plot(abs(rhs_z), w_vec,'o');
xlabel('absolute  res')
ylabel('weight')

end

 function [x_sol, rhs_x, rhs_z, Ax, Az, A,rhs,Apw,rhs_w,wHub_z,R]= getGPSgraph(x_vec,Y,n,x0,sig_e,Sp_sq,n_sat,n_state,x_sat,y_sat,z_sat,c,tm)
% x_vec: linearization point for the state vector
% x0: vehicle states at time 0

%state vector into n_state by n matrix
x1=reshape(x_vec,[n_state,n]); 
% rhs vector for x
for i=1:n
    G=zeros(n_state,n_state);
    if i==1
        T=tm(1,i);
        F2=[1 T; 0 1];
        G(4:n_state,4:n_state)=F2;
        a_vec_i=x1(:,i)-G*x0;
        rhs_x=Sp_sq*a_vec_i;
    else
        T=tm(1,i)-tm(1,i-1);
        F2=[1 T; 0 1];
        G(4:n_state,4:n_state)=F2;      
        a_vec_i=x1(:,i)-G*x1(:,i-1);
        rhs_x=[rhs_x; Sp_sq*a_vec_i];
    end
end
% rhs vector for z
rhs_z=[]; 
for i=1:n
  v_states=x1(:,i);
  r_meas=Y(:,i);
  s_coord=[x_sat(:,i) y_sat(:,i) z_sat(:,i)];
  [c_vec_i,~]=getRhs_z(s_coord,v_states,r_meas,c,n_sat(i));  
  Se_sq=(inv(eye(n_sat(i))*sig_e))^0.5;
  rhs_z=[rhs_z; Se_sq*c_vec_i];  
end
% combine rhs's. % truncate rhs_x.  omit first state observation (start with n_state+1)
rhs=[rhs_x(n_state+1:end); rhs_z];
% A matrix for x
% diagonals
Ax1=cell(n,1);
for i=1:n
    Ax1{i}= sparse(-Sp_sq*eye(n_state));
end
Ax=blkdiag(Ax1{:});
%lower diagonals
for i=2:n
    G=zeros(n_state,n_state);
    T=tm(1,i)-tm(1,i-1);
    F2=[1 T; 0 1];
    G(4:n_state,4:n_state)=F2;
    Ax(n_state*(i-1)+[1:n_state],n_state*(i-2)+[1:n_state])=Sp_sq*G; 
end
% A matrix for z
Az1=cell(n,1);    
for i=1:n
  v_states=x1(:,i);
  s_coord=[x_sat(:,i) y_sat(:,i) z_sat(:,i)];
  Hi=getH(s_coord,v_states,c,n_sat(i));
  Se_sq=(inv(eye(n_sat(i))*sig_e))^0.5;
  Az1{i}= sparse(Se_sq*Hi);
end
Az=blkdiag(Az1{:});
% combine the A matrices (start Ax with n_state+1))
A=[Ax(n_state+1:end,:); Az];
%% factor graph approach
% robust estimate of scale- (Montgomery et al., 2006, p. 374)
% apply to residuals of measurements only
s_robust=median(abs(rhs_z-median(rhs_z)))/0.6745;
z_res=rhs_z/s_robust;
% Ramsay's Ea function with a=0.3
%Ra=0.3;
%wHub_z=exp(-Ra*abs(z_res'));
% Welsch function 
%Rw=2;
%wHub_z=exp(-(z_res'/Rw).^2);
% Geman-McClure weight function (Agarwal, 2013)
%k_gm=4;
%wHub_z=1./(1+(z_res'/k_gm).^2).^2;
% Cauchy function
k_c=3.5;
wHub_z=1./(1+(z_res'/k_c).^2);
% Huber's t with t=2
%t_hub=1.345;
% for i=1:length(z_res)
%     if abs(z_res(i))<=t_hub
%         wHub_z(i)=1;
%     else
%         wHub_z(i)=t_hub/abs(z_res(i));
%     end
% end
% % Andrews' wave with a=1.339
% for i=1:length(z_res)
%     if abs(z_res(i))<=1.339*pi
%         wHub_z(i)=sin(z_res(i)/1.339)/(z_res(i)/1.339);
%     else
%         wHub_z(i)=0;
%     end
% end
% % Hampbels' 17A
% for i=1:length(z_res)
%     if abs(z_res(i))    <=1.7
%         wHub_z(i)=1;
%     elseif abs(z_res(i))<=3.4
%         wHub_z(i)=1.7/abs(z_res(i));
%     elseif abs(z_res(i))<=8.5
%         wHub_z(i)=(1.7*(8.5-abs(z_res(i))))/(abs(z_res(i))*(8.5-3.4));
%     else
%         wHub_z(i)=0;
%     end
% end
weights=[ones(n_state*(n-1),1);wHub_z'];
nw=length(weights);

%%% apply to residuals of measurements and process - Ramsay's Ea function
% s_robust=median(abs(rhs-median(rhs)))/0.6745;
% z_res=rhs/s_robust;
% Ra=0.25;
% wHub_z=exp(-Ra*abs(z_res'));
% weights=wHub_z';
%%%

% diagonal weights matrix
Wsq=sparse([1:nw],[1:nw],sqrt(weights),nw,nw);
% find weighted A matrix and rhs vector 
Aw=Wsq*A;
rhs_w=Wsq*rhs;
% reorder columns of A to exploit sparsity
indcolA=colamd(Aw);
Apw=Aw(:,indcolA);
% inverse reordering to find original variables
[~,indInvSortA]=sort(indcolA);
%QR factorization 
%xp=Apw\rhs_w; 
Rdum = qr(Apw,0);
R=triu(Rdum); 
xp = R\(R'\(Apw'*rhs_w));
% reorder rows of solution vector 
x_sol_factor=xp(indInvSortA);
x_sol=x_sol_factor;
end

function [c_vec,hvec]=getRhs_z(s_coord,v_states,r_meas,c,nsat_i) 
% s_coord:  n_sat by 3 satellite positions
% v_states: 1 by 5 vehicle state estimates ,  x y z b(clock error), d(clock drift) 
% r_meas: n_sat by 1, pseudorange 
% h(i)=r_i+c*b, noiseless pseudorange to satellite i, where c speed of sound and r_i=sqrt((sat_xi-x)^2+()^2+()^2)
% for h(i) see Bar Shalom, p. 506

for i=1:nsat_i
    hvec(i,1)=sqrt((s_coord(i,1)-v_states(1))^2+...
        (s_coord(i,2)-v_states(2))^2+(s_coord(i,3)-v_states(3))^2)+c*v_states(4)+(7.3e-5/3e8)*(s_coord(i,1)*v_states(2)-s_coord(i,2)*v_states(1));
end
c_vec=r_meas(1:nsat_i)-hvec;
end

function H=getH(s_coord,v_states,c,nsat_i)
n_state=length(v_states);
H=zeros(nsat_i,n_state);
% see Bar Shalom, p. 513
for i=1:nsat_i
    hi=sqrt((s_coord(i,1)-v_states(1))^2+(s_coord(i,2)-v_states(2))^2+(s_coord(i,3)-v_states(3))^2);
    Hi1=-1/hi*[s_coord(i,1)-v_states(1) s_coord(i,2)-v_states(2) s_coord(i,3)-v_states(3)]+...
        [-(7.3e-5/3e8)*s_coord(i,2) (7.3e-5/3e8)*s_coord(i,1) 0];
    H(i,:)=[Hi1 c 0];
end
end

function [sig_pHat1_std, sig_pHat2_std, sig_eHat_std]=getUnbVarEstSLAM(rx1,rx2,rz,Ax1,Ax2,Az,flag_scale_cov_change)

A1=Ax1;
A2=Ax2;
A3=Az;
r1=rx1;
r2=rx2;
r3=rz;
mx1=size(A1,1);
mx2=size(A2,1);
mz=size(A3,1);
p=size(A1,2);
%QR decompositon of A
Rdum = qr([A1;A2;A3],0);    % full matrix
R=triu(Rdum);               % upper triangular entries   

% find traces of H1, H2, H3, H1^2, etc. without matrix inversion
M1=A1/R;
clear A1
tH1=sum(diag(M1*M1'));
tH11=sum(sumsqr(M1*M1'));
M2=A2/R;
clear A2
tH2=sum(diag(M2*M2'));
tH22=sum(sumsqr(M2*M2'));
tH1221=sum(diag((M1'*M1)*(M2'*M2)));
clear M2 M1
%M3=A3/R;
%clear A3
%tH3=sum(diag(M3*M3'));
%tH33=sum(sumsqr(M3*M3'));
%clear M3

% usually A3 has the largest number of rows than A1 and A2. Find traces of
% A3 by subtraction
tH3=p-tH1-tH2;

tD11_1=mx1-2*tH1+tH11;
%tD22_1=0.5*(tH1+tH2-tH3-(tH11+tH22-tH33));
tD22_1=tH1221;

tD11_2=tD22_1;
tD22_2=mx2-2*tH2+tH22;
%tD33_2=-0.5*(tH1-tH2-tH3-(tH11-tH22-tH33));
tH2332=tH2-tH22-tH1221;
tD33_2=tH2332;
tD22_3=tD33_2;

%tD33_1=0.5*(tH1-tH2+tH3-(tH11-tH22+tH33));
tH1331=tH1-tH11-tH1221;
tD33_1=tH1331;
tD11_3=tD33_1;
%tD33_3=mz-2*tH3+tH33;
tD33_3=mz-2*tH3+(tH3-tH1331-tH2332);

% estimate three variance components
Coef_mat=[  tD11_1 tD22_1 tD33_1;
            tD11_2 tD22_2 tD33_2;
            tD11_3 tD22_3 tD33_3];

%var_est=Coef_mat\[r1'*r1; r2'*r2; r3'*r3];

if flag_scale_cov_change
    true_res = [r1.' * r1; r2.' * r2; r3.' * r3];
    %Compute (the log of) what residuals should be if it had converged to 
    %the right covariances
    des_res = log(Coef_mat * ones(3,1));
    curr_res = true_res;
    var_est = Coef_mat\curr_res;
    %Loop until we don't move too far
    %Too far is arbitrarily order of magnitude...
    while nnz(var_est > 5.0) > 0 || nnz(var_est < 0.2) > 0
        %Move "curr_res" towards the des_res if moved too far
        curr_res = exp((des_res + log(curr_res))/2);
        var_est = Coef_mat\curr_res;    
    end
else
    var_est=Coef_mat\[r1'*r1; r2'*r2; r3'*r3];
end


sig_pHat1_std=var_est(1);
sig_pHat2_std=var_est(2);
sig_eHat_std=var_est(3);

end

