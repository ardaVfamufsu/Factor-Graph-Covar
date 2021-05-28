% Object moving on a plane with constant velocity. No landmarks. Measurements are GPS x,y coordinates 
% Process noise has 2 variance components, sig_p1 for x coordinates and sig_p2 for y coordinates
% Q=[sig_p1 0 0 0; 0 sig_p1 0 0; 0 0 sig_p2 0; 0 0 0 sig_p2];
% Measurement variance sig_e is the same for both x and y coordinates.
% ==== user inputs - comment some if running in batch mode
close all;
clear;
% flag =1 for unbiased, 0 for biased variance estimators
flag_unbiased=1;
flag_new=1;  %Should hopefully be correct this time...
% flag =1 for robust (M), 0 for nonrobust state estimators (biased or unbiased variance estimators)
flag_robust=0;
% flag =1 for errors with outliers, 0 for no outliers
flag_outliers=0;
% flag=1 if trying to "scale" adjusted covariance to keep with order of magnitude change
flag_scale_cov_change=1;
% number of MC runs
n_MC=100;
% number of times the robot will repeat the cycles
n_rep_cycle=1;
% random number seed
% rng(0);
% ==== end of user inputs
% sampling period
T=1;
% known initial vehicle state: (x xd y yd)
x0=[0 2  0  0]';
% process dynamics and measurement matrices
G=[1 T 0 0; 0 1 0 0; 0 0 1 T;0 0 0 1];
F=[1 0 0 0; 0 0 1 0];
% controls on velocities - S shape motion
vel_x=2; vel_y=2;
control_x_vel_cycle=[0 0 0 0 -vel_x 0 0 0 0  vel_x 0 0 0 0 -vel_x 0 0 0 0 vel_x];
control_y_vel_cycle=[0 0 0 0  vel_y 0 0 0 0 -vel_y 0 0 0 0 -vel_y 0 0 0 0 vel_y];
%length of each cycle
n_cycle=length(control_x_vel_cycle);
control_mat_cycle=zeros(4,n_cycle);
control_mat_cycle(2,:)=control_x_vel_cycle;
control_mat_cycle(4,:)=control_y_vel_cycle;
% replicated cycles
control_mat=repmat(control_mat_cycle,1,n_rep_cycle);
n=size(control_mat,2);
% number of states, number of measurements
n_state=4;
n_meas=2;
% no of rows and columns (no of parameters) in A matrix
m=(n_state+n_meas)*n;
p=n_state*n;
% no of rows for process and measurement blocks
mx=p;
mz=m-mx;
% indexes of states that have variance sig_p1 and sig_p2, respectively
ind_x1=union([1:n_state:mx],[2:n_state:mx]);
ind_x2=union([3:n_state:mx],[4:n_state:mx]);
mx1=length(ind_x1);
mx2=length(ind_x2);
% measurement and process noise, variance
sig_e=1.5;
sig_p1=0.5;     % x coordinate and x velocity
sig_p2=0.2;     % y coordinate and y velocity. 
% measurement error mixture parameters (for outliers)
w1=0.90; w2=1-w1;
mu_e1=0;mu_e2=0;
sig_e1=sig_e; sig_e2=10;
% IRLS parameters
max_err=0.001;
if flag_unbiased
    max_iter=50;    % use with unbiased estimates
else
    max_iter=3;     % use with biased estimates
end
% modify for single realization plots, no outliers
if n_MC==1 && flag_outliers==0
    sig_e=0.5;sig_p1=0.05;sig_p2=0.03;
    max_err=0.000001;
    rng(4);
end
% modify for single realization plots, with outliers
if n_MC==1 && flag_outliers==1
    sig_e=0.5;sig_p1=0.05;sig_p2=0.03;
    sig_e1=sig_e; sig_e2=5;
    max_err=0.000001;
    rng(4);
end
% measurement and process noise, variance matrices
V=sig_e*eye(n_meas);
W=[sig_p1 0 0 0; 0 sig_p1 0 0; 0 0 sig_p2 0; 0 0 0 sig_p2];
% intialized values of measure and process noise, for IRLS
sig_eHat0=1; sig_pHat0=1;
%% state estimation, linearization approach - See Dellaert and Kaess  (2006)
% number of linearizations
m_iter=5;
% mahalanobis distance storage vector
distMah_irls=zeros(n_MC,1); 
distMah_irls2=zeros(n_MC,1); 
for j=1:n_MC
    % process noise
    w=mvnrnd(zeros(n,4),W,n);
    if flag_outliers 
        % measure noise: mixture components, x and y coordinates sampled together from either component
        flag_m1=rand(n,1)<w1; 
        flag_m2=1-flag_m1;
        flag_m11=kron(flag_m1,ones(n_meas,1));
        flag_m22=kron(flag_m2,ones(n_meas,1));
        ind_m1=find(flag_m11); 
        ind_m2=find(flag_m22);
        e_v=zeros(n_meas*n,1);
        e_v(ind_m1)=sig_e1*randn(length(ind_m1),1)+mu_e1;
        e_v(ind_m2)=sig_e2*randn(length(ind_m2),1)+mu_e2;
        e=transpose(reshape(e_v,n_meas,n));  
    else
        % measure noise: no outliers
        e=mvnrnd(zeros(n,n_meas),V,n);
    end
    % find states from the process model
    x=G*x0+control_mat(:,1)+transpose(w(1,:));
    for i=2:n
        x(:,i)=G*x(:,i-1)+control_mat(:,i)+transpose(w(i,:));
    end
    x1_true=x(1,:);
    x2_true=x(3,:);
    % measurements
    Yn(:,1)=x1_true'+e(:,1);      
    Yn(:,2)=x2_true'+e(:,2); 
    % IRLS iteration parameters
    err_sigSol=99;
    c1=1;
    % linearization point
    x_vec0=[repmat(x0,n,1)];
    x_sol_prev=x_vec0;
    x_sol=x_vec0;    
    % initialize process and measurement variances
    sig_eHat=sig_eHat0; sig_p1Hat=sig_pHat0; sig_p2Hat=sig_pHat0;
    % initialize condition number
    cond_noA=1e-6;      
    % cputime start
    cptime0=cputime;
    % iteratively reweighted LS (Montgomery et al., 2006, p. 374)
    while c1<max_iter &&err_sigSol>max_err&&cond_noA<100
        %measurement and process covariances and inverse square roots
        V_est=sig_eHat*eye(n_meas);
        W_est=[sig_p1Hat 0 0 0; 0 sig_p1Hat 0 0; 0 0 sig_p2Hat 0; 0 0 0 sig_p2Hat];
        Se_sq_est=(inv(V_est))^0.5;
        Sp_sq_est=(inv(W_est))^0.5;
        for i=1:m_iter           
            if flag_robust
                [del_x_sol,~, ~,~, ~,Aunw,~,Aw,rhs_w,~,~,~,R,Apw]=linearSLAM_noLandmark_gps_Mestimate(x_vec0,Yn',control_mat,n,G,F,x0,Se_sq_est,Sp_sq_est); 
                A=Apw; % weighted and permuted A matrix
                Aunpermute=Aw; % weighted A matrix
                rhs_x=rhs_w(1:mx);
                rhs_z=rhs_w(mx+1:m);
            else
                if flag_new
                    [del_x_sol,rhs_x,rhs_z,A,unw_A,W,unw_R] = linearSLAM_newWeighting(x_vec0,Yn',control_mat,n,G,F,x0,Se_sq_est,Sp_sq_est);
           
                else
                    % A and Aunpermuted: A and unpermuted A matrices
                    [del_x_sol,rhs_x, rhs_z,~, ~,Aunpermute,~,R,A]=linearSLAM_noLandmark_gps(x_vec0,Yn',control_mat,n,G,F,x0,Se_sq_est,Sp_sq_est);   
                end
            end
            x_vec0=x_vec0+del_x_sol;  
            % linearization error
            lin_err(1,i)=dot(del_x_sol,del_x_sol);
        end
        x_sol=x_vec0;
        x_sol_prev=x_sol;
        rhs_x1=rhs_x(ind_x1);
        rhs_x2=rhs_x(ind_x2);
        Ax1=A(ind_x1,:);
        Ax2=A(ind_x2,:);
        Az=A(mx+1:m,:); 
        if flag_unbiased
            if flag_new
                idx{1}=ind_x1;
                idx{2}=ind_x2;
                idx{3}=mx+1:m;
                cov_est = getUnbVarEstSLAM_new(unw_A,W,unw_R,idx);
                sig_p1Hat_new = cov_est{1};
                sig_p2Hat_new = cov_est{2};
                sig_eHat_new = cov_est{3};
                
            else
                % unbiased estimate of standardized measurement, process and combined variances
                [sig_p1Hat_std,sig_p2Hat_std, sig_eHat_std]=...
                    getUnbVarEstSLAM_fastQR_3Components(R,rhs_x1,rhs_x2,rhs_z,Ax1,Ax2,Az,flag_scale_cov_change); % using Fast QR
            end
        else
            % biased (nonrobust) estimate of standardized measurement, process and combined variances
            sig_p1Hat_std=1/(mx1-1)*dot(rhs_x1,rhs_x1);  
            sig_p2Hat_std=1/(mx2-1)*dot(rhs_x2,rhs_x2);  
            sig_eHat_std=1/(mz-1)*dot(rhs_z,rhs_z);
        end
        % update previous values
        sig_p1Hat_prev=sig_p1Hat;
        sig_p2Hat_prev=sig_p2Hat;
        sig_eHat_prev=sig_eHat;  
        if flag_new
           sig_p1Hat =sig_p1Hat_new;
           sig_p2Hat =sig_p2Hat_new;
           sig_eHat =sig_eHat_new;
        else
            % find new (unstandardized) variances -using unbiased estimates 
            sig_p1Hat=(sig_p1Hat_std)*sig_p1Hat_prev;     
            sig_p2Hat=(sig_p2Hat_std)*sig_p2Hat_prev;  
            sig_eHat=(sig_eHat_std)*sig_eHat_prev;
        end
        % error vector for variance estimate iterations
        err_vec=[sig_eHat-sig_eHat_prev; sig_p1Hat-sig_p1Hat_prev; sig_p2Hat-sig_p2Hat_prev];
        err_sigSol=dot(err_vec,err_vec);            
        % iteration counter
        c1=c1+1;      
        S_st(c1,:)=[sig_eHat sig_p1Hat sig_p2Hat];        
        % condition number of weighted, permuted A matrix
        cond_noA=cond(A);   
        % for use in n_MC=1
        if n_MC==1
            disp(['iteration ',num2str(c1)]);
            disp('process 1, process 2 and measure variances')
            disp([sig_p1Hat,sig_p2Hat, sig_eHat])
            disp('process 1, process 2 and measure standardized variances')
            disp([sig_p1Hat_std, sig_p2Hat_std,sig_eHat_std])   
            disp(['linearization state sum of squared errors in ', num2str(m_iter),' iterations'])
            disp(lin_err) 
        end
    end
    % cputime for state and covariance estimation
    cptime(j)=cputime-cptime0;
    %j
    
    %% At this point, modify to make the variance estimates be conservative
%     dof = 3
    % store estimates
    var_vect(:,j)=[sig_eHat;sig_p1Hat;sig_p2Hat];
    % error vector from true trajectory
    v_est_irls=reshape(x_sol,[4,n]);
    e_irls=[v_est_irls(1,:)-x1_true; v_est_irls(3,:)-x2_true]; 
    % euclidian distance from true trajectory
    distSq_irls(j)=dot(e_irls(:),e_irls(:));    
    % time trace of state covariances
    
    %Compute Mahalanobis distance for e_irls
    Cx=inv(Aunpermute'*Aunpermute); % estimated covar
    C_pos = Cx(1:2:end,1:2:end);
    distMah_irls2(j) = e_irls(:).'*inv(C_pos)*e_irls(:);
    for t=1:n
         % extract [11 13; 31 33] portion of each 4 by 4 block diagonal term of Cx
         r1=n_state*(t-1)+1;
         Cy_vec(:,:,t)=[Cx(r1,r1) Cx(r1,r1+2);Cx(r1,r1+2) Cx(r1+2,r1+2)];
         % mahalanobis distance from true trajectory
         distMah_irls(j)=distMah_irls(j)+transpose(e_irls(:,t))*inv(Cy_vec(:,:,t))*e_irls(:,t);
    end 
    if distMah_irls(j) > 500
        plot(e_irls)
    end
end

vest_mean=mean(var_vect');
var_true=[sig_e;sig_p1;sig_p2];
bias_estimates=dot(vest_mean'-var_true,vest_mean'-var_true);
disp('mean of the variance estimates');
disp(vest_mean);
disp('true variances ');
disp(var_true');
disp('bias of the estimates ')
disp(bias_estimates)
disp('average cputime (sec)')
disp(mean(cptime))
%% output results of estimation (when n_MC=1)
if n_MC==1
%     % convergence of estimates
%     figure;
%     plot([1:c1-1],S_st(2:end,1),'d-k'); hold on; plot([1 c1-1],[sig_e sig_e],'r');xlabel('iter');ylabel('Measure var, \sigma^2_r');
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     %saveas(gcf,'varConverge1.jpg')
%     figure;
%     plot([1:c1-1],S_st(2:end,2),'d-k');hold on; plot([1 c1-1],[sig_p1 sig_p1],'r');xlabel('iter');ylabel('Process var 1, \sigma^2_{q1}');
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     figure;
%     plot([1:c1-1],S_st(2:end,3),'d-k');hold on; plot([1 c1-1],[sig_p2 sig_p2],'r');xlabel('iter');ylabel('Process var 2, \sigma^2_{q2}');
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     %saveas(gcf,'varConverge2.jpg')
    % vehicle trajectory and measurements
    figure;    
    y_noNoise=[zeros(1,5) 2:2:10 10*ones(1,5) 8:-2:0];
    x_noNoise=[2:2:10 10*ones(1,5) 12:2:20 20*ones(1,5)];
    plot(x_noNoise,y_noNoise,'k-^','Linewidth',1);
    hold on;
    plot(x1_true,x2_true,'o-b','Linewidth',1);
    xlabel('x loc')
    ylabel('y loc')
    plot(Yn(:,1),Yn(:,2),'d-r','Linewidth',1);
    legend('control','true trajectory','measured');
    axis([0 30 -2 20]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gca,'FontSize',16)
    saveas(gcf,'measSim.jpg')
    %vehicle trajectory and estimates
    v_est=reshape(x_sol,[4,n]);
    figure;
    plot(v_est(1,:),v_est(3,:),'dg-','Linewidth',1);
    hold on;
    plot(x1_true,x2_true,'o-b','Linewidth',1);
    xlabel('x loc')
    ylabel('y loc')
    legend('Estimated Trajectory','True Trajectory');
    axis([0 30 -2 20]);
    fig1=gcf;
    % time trace of state covariances- plot conf ellipses
    Cx=inv(Aunpermute'*Aunpermute); % estimated covar
    for i=1:n
        % extract [11 13; 31 33] portion of each 4 by 4 block diagonal term of Cx
        r1=4*(i-1)+1;
        Cmat=[Cx(r1,r1) Cx(r1,r1+2);Cx(r1,r1+2) Cx(r1+2,r1+2)];
        Cy_vec(:,:,i)=Cmat;
    end
    % plot conf ellipses on trajectory
    f2=copyobj(fig1,0);
    for i=1:n
        [Xellipse,Yellipse,eigvec,eigval]=findConfEllipse(Cy_vec(:,:,i),[v_est(1,i); v_est(3,i)]);
        plot(Xellipse,Yellipse,'g','Linewidth',1,'HandleVisibility','off'); 
        quiver(v_est(1,i), v_est(3,i),3*sqrt(eigval(1,1))*eigvec(1,1),3*sqrt(eigval(1,1))*eigvec(2,1),'g','HandleVisibility','off')
    end
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gca,'FontSize',16);
    saveas(gcf,'estState.jpg')
    disp(['Mah Dist is ',num2str(distMah_irls)]);
end

%% output monte carlo results (when n_MC>1)
if n_MC>1
    e_above = nnz(var_vect(1,:)>sig_e)
    e_below = length(var_vect)-e_above;
    p1_above = nnz(var_vect(2,:)>sig_p1)
    p2_above = nnz(var_vect(3,:)>sig_p2)
    
    
    % histograms of error variance estimates and true values
    figure;
    histogram(var_vect(1,:),[-1.25:0.2:4],'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_e sig_e],[0 300],'b--','LineWidth',2);xlabel('Measure Var, \sigma^2_r');ylabel('Frequency')
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'varHistR.jpg')
    %title(['measure var estimates, 2.5, 50, 97.5%: ',num2str(round(prctile(var_vect(1,:),[2.5 50 97.5]),3))])
    figure;
    histogram(var_vect(2,:),'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_p1 sig_p1],[0 n_MC/4],'b--','LineWidth',2);xlabel('Process Var,\sigma^2_{q1}');ylabel('Frequency')
    %title(['proc var estimates, 2.5, 50, 97.5%: ',num2str(round(prctile(var_vect(2,:),[2.5 50 97.5]),3))])
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));    
    saveas(gcf,'varHistQ1.jpg')
    figure;
    histogram(var_vect(3,:),'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_p2 sig_p2],[0 400],'b--','LineWidth',2);xlabel('Process Var,\sigma^2_{q2}');ylabel('Frequency')
    %title(['proc var estimates, 2.5, 50, 97.5%: ',num2str(round(prctile(var_vect(3,:),[2.5 50 97.5]),3))])
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'varHistQ2.jpg')
    % histograms of position estimate error (Mahalanobis)
    figure;
    histogram(distMah_irls,[0:20:360],'FaceColor','w','EdgeColor','k');hold on;
    xlabel('Mahalanobis Distance of State Est Error');ylabel('Frequency');
    %title('Unbiased Estimate of Variances')
    %plot(prctile(distMah_irls,2.5)*ones(1,2),[0 300],'g--');
    %plot(prctile(distMah_irls,50)*ones(1,2),[0 300],'g--');
    %plot(prctile(distMah_irls,97.5)*ones(1,2),[0 300],'g--');
    plot([2*n 2*n],[0 300],'--r','LineWidth',2)
    axis([0 400 0 300]);
    set(gca,'FontSize',16);
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'MahHist.jpg');
    disp(['Mah Dist 2.5, 50, 97.5 Percentiles: ',num2str(round(prctile(distMah_irls,[2.5 50 97.5]),3))])

    % Mah try 2
    figure;
    histogram(distMah_irls2,[0:20:360],'FaceColor','w','EdgeColor','k');hold on;
    xlabel('Mahalanobis Distance of State Est Error2');ylabel('Frequency');
    %title('Unbiased Estimate of Variances')
    %plot(prctile(distMah_irls,2.5)*ones(1,2),[0 300],'g--');
    %plot(prctile(distMah_irls,50)*ones(1,2),[0 300],'g--');
    %plot(prctile(distMah_irls,97.5)*ones(1,2),[0 300],'g--');
    plot([2*n 2*n],[0 300],'--r','LineWidth',2)
    axis([0 400 0 300]);
    set(gca,'FontSize',16);
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'MahHist.jpg');
    disp(['Mah Dist 2.5, 50, 97.5 Percentiles: ',num2str(round(prctile(distMah_irls2,[2.5 50 97.5]),3))])

end