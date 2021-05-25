% forward velocity and bearing model- nonlinear process equation
% ==== user inputs - comment some if running in batch mode
clear;
close all;
% number of times the robot will repeat the cycles
n_rep_cycle=1;
% robust estimation (robust implementation not yet ready)
flag_robust=0;
% unbiased estimation
flag_unbiased=0;
% flag=1 if trying to "scale" adjusted covariance to keep with order of magnitude change
flag_scale_cov_change=1;
% number of MC runs
n_MC=100;
% random number seed
rng(2);
% ==== end of user inputs
% sampling period
T=1;
% known initial vehicle state: (x y theta v)
x0=[0 0 0 2]';
% measurement matrices
F=[1 0 0 0; 0 1 0 0];
% controls on velocities - forward and rotation
vel=1; 
vel_t=pi/2;
% more nonlinear patterns of angular velocity control
control_t_vel_cycle=[0 0 0 0  vel_t 0 0 0 0 -vel_t 0 0 0 0 -vel_t 0 0 0 0 vel_t];
%control_t_vel_cycle=[0 0 vel_t 0  0 -vel_t 0 0 -vel_t 0 0 vel_t 0 0 vel_t 0  0 -vel_t 0 0 -vel_t 0 0 vel_t];
%control_t_vel_cycle=[0 vel_t 0  -vel_t 0 -vel_t 0 vel_t 0 vel_t 0 -vel_t 0 -vel_t 0 vel_t 0 vel_t 0  -vel_t 0 -vel_t 0 vel_t 0 vel_t 0 -vel_t 0 -vel_t 0 vel_t];

%length of each cycle
n_cycle=length(control_t_vel_cycle);
control_vel_cycle=zeros(1,n_cycle);
control_mat_cycle=zeros(4,n_cycle);
control_mat_cycle(4,:)=control_vel_cycle;
control_mat_cycle(3,:)=control_t_vel_cycle;
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
% measurement and process noise, variance
sig_e=1.5;
sig_p1=0.5;     
sig_p2=0.2;
% indexes of states that have variance sig_p1 and sig_p2, respectively
ind_x1=union([1:n_state:mx],[2:n_state:mx]);
ind_x2=union([3:n_state:mx],[4:n_state:mx]);
mx1=length(ind_x1);
mx2=length(ind_x2);
% number of linearizations
m_iter=5;
% IRLS parameters
max_err=0.001;
if flag_unbiased
    max_iter=50;    % use with unbiased estimates
else
    max_iter=3;     % use with biased estimates
end
% modify for single realization plots, no outliers
if n_MC==1
    sig_e=0.5;sig_p1=0.05;sig_p2=0.03;
    max_err=0.000001;
    rng(4);
end
% measurement and process noise, variance matrices
V=sig_e*eye(n_meas);
W=[sig_p1 0 0 0; 0 sig_p1 0 0; 0 0 sig_p2 0; 0 0 0 sig_p2];
%% MC simulations
% mahalanobis distance storage vector
distMah_irls=zeros(n_MC,1); 
ind_large=[];
for j=1:n_MC
    % process noise
    w=mvnrnd(zeros(n,n_state),W,n);
    % measure noise: no outliers
    e=mvnrnd(zeros(n,n_meas),V,n);
    % find states from the process model
    for i=1:n
        if i==1
            x_prev=x0(1);
            y_prev=x0(2);
            theta_prev=x0(3);
            v_prev=x0(4);
        else
            x_prev=x(1,i-1);
            y_prev=x(2,i-1);
            theta_prev=x(3,i-1);
            v_prev=x(4,i-1);        
        end
        x(:,i)=[x_prev;y_prev;theta_prev;v_prev]+...
        [v_prev*T*cos(theta_prev);v_prev*T*sin(theta_prev); 0; 0]+control_mat(:,i)+transpose(w(i,:));
    end
    x1_true=x(1,:);
    x2_true=x(2,:);
    % measurements
    Yn(:,1)=x1_true'+e(:,1);      
    Yn(:,2)=x2_true'+e(:,2); 
    % IRLS iteration parameters
    err_sigSol=99;
    c1=1;
    sig_eHat0=1; sig_pHat0=1;
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
                [del_x_sol,~, ~,~, ~,Aunw,~,Aw,rhs_w,~,~,~,R,Apw]=TEST_forwVel_linearSLAM_noLandmark_gps_Mestimate(x_vec0,Yn',control_mat,n,F,x0,Se_sq_est,Sp_sq_est,T);
                A=Apw; % weighted and permuted A matrix
                Aunpermute=Aw; % weighted A matrix
                rhs_x=rhs_w(1:mx);
                rhs_z=rhs_w(mx+1:m);
            else               
                % A and Aunpermuted: A and unpermuted A matrices
                [del_x_sol,rhs_x, rhs_z,~, ~,Aunpermute,rhs,R,A]=linearSLAM_forwardVel_noLandmark_gps(x_vec0,Yn',control_mat,n,F,x0,Se_sq_est,Sp_sq_est,T);   
            end
            %x_vec0=x_vec0+del_x_sol;                          
            %check linearity
            scalar = 1;
            test_ratio = 0;	
            while (scalar > 1E-5 & ~(test_ratio<4 & test_ratio > .25))
                prop_vec = x_vec0+del_x_sol*scalar;
                %Check the cost to make sure the linearization assumption was valid
                curr_cost = dot(rhs, rhs);
                pred_z = rhs - A * del_x_sol*scalar;
                pred_cost = dot(pred_z, pred_z);
                [~,~, ~,~, ~,~,new_z,~,~]=linearSLAM_forwardVel_noLandmark_gps(prop_vec,Yn',control_mat,n,F,x0,Se_sq_est,Sp_sq_est,T);   
                true_cost = dot(new_z, new_z);
                test_ratio = abs(curr_cost - pred_cost)/abs(curr_cost - true_cost);
                scalar = scalar / 2.0;
            end
            x_vec0 = prop_vec;
            %if scalar < 1E-5
            %    warning('Scaling the Gauss Newton more than desired.  Are your Jacobians correct?');
            %end
            %    fprintf ('Iteration %d, scaled by %f\n',i,scalar * 2.0)           
            %check linearity- end
                       
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
            % unbiased estimate of standardized measurement, process and combined variances
            [sig_p1Hat_std,sig_p2Hat_std, sig_eHat_std]=...
                getUnbVarEstSLAM_fastQR_3Components(R,rhs_x1,rhs_x2,rhs_z,Ax1,Ax2,Az,flag_scale_cov_change); % using Fast QR  
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
        % find new (unstandardized) variances -using unbiased estimates 
        sig_p1Hat=(sig_p1Hat_std)*sig_p1Hat_prev;     
        sig_p2Hat=(sig_p2Hat_std)*sig_p2Hat_prev;  
        sig_eHat=(sig_eHat_std)*sig_eHat_prev;
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
    % store estimates
    var_vect(:,j)=[sig_eHat;sig_p1Hat;sig_p2Hat];
    % index of too large values
    if max(abs(var_vect(:,j)))>100 ind_large=[ind_large j]; end;
    % error vector from true trajectory
    v_est_irls=reshape(x_sol,[4,n]);
    e_irls=[v_est_irls(1,:)-x1_true; v_est_irls(2,:)-x2_true]; 
    % euclidian distance from true trajectory
    distSq_irls(j)=dot(e_irls(:),e_irls(:));    
    % time trace of state covariances
    Cx=inv(Aunpermute'*Aunpermute); % estimated covar
    for t=1:n
         % extract [11 12; 21 22] portion of each 4 by 4 block diagonal term of Cx
         r1=n_state*(t-1)+1;
         Cy_vec(:,:,t)=[Cx(r1,r1) Cx(r1,r1+1);Cx(r1,r1+1) Cx(r1+1,r1+1)];
         % mahalanobis distance from true trajectory
         distMah_irls(j)=distMah_irls(j)+transpose(e_irls(:,t))*inv(Cy_vec(:,:,t))*e_irls(:,t);
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
    plot(v_est(1,:),v_est(2,:),'dg-','Linewidth',1);
    hold on;
    plot(x1_true,x2_true,'o-b','Linewidth',1);
    xlabel('x loc')
    ylabel('y loc')
    legend('Estimated Trajectory','True Trajectory');
    
    fig1=gcf;
    % time trace of state covariances- plot conf ellipses
    Cx=inv(Aunpermute'*Aunpermute); % estimated covar
    for i=1:n
        % extract [11 12; 21 22] portion of each 4 by 4 block diagonal term of Cx
        r1=n_state*(i-1)+1;
        Cmat=[Cx(r1,r1) Cx(r1,r1+1);Cx(r1,r1+1) Cx(r1+1,r1+1)];
        Cy_vec(:,:,i)=Cmat;
    end
    % plot conf ellipses on trajectory
    f2=copyobj(fig1,0);
    hold on;
    for i=1:n
        [Xellipse,Yellipse,eigvec,eigval]=findConfEllipse(Cy_vec(:,:,i),[v_est(1,i); v_est(2,i)]);
        plot(Xellipse,Yellipse,'g','Linewidth',1,'HandleVisibility','off'); 
        quiver(v_est(1,i), v_est(2,i),3*sqrt(eigval(1,1))*eigvec(1,1),3*sqrt(eigval(1,1))*eigvec(2,1),'g','HandleVisibility','off')
    end
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gca,'FontSize',16);
    %saveas(gcf,'estState.jpg')
    disp(['Mah Dist is ',num2str(distMah_irls)]);
    figure;
    plot(v_est(3,:));
    hold on;
    stairs(x0(3)+cumsum(control_mat(3,:)),'m')
    legend('Estimated Rotation speed','True ');
    
    figure;
    plot(v_est(4,:));
    hold on;
    stairs(x0(4)+cumsum(control_mat(4,:)),'m')
    legend('Estimated Forward speed','True ');
end
%% output monte carlo results (when n_MC>1)
if n_MC>1
    % omit large variance values
    var_vect(:,ind_large)=[];
    distMah_irls(ind_large)=[];
    % histograms of error variance estimates and true values
    figure;
    histogram(var_vect(1,:),[-1.25:0.2:4],'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_e sig_e],[0 300],'b--','LineWidth',2);xlabel('Measure Var, \sigma^2_r');ylabel('Frequency')
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'varHistRNonlinear.jpg')
    figure;
    histogram(var_vect(2,:),'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_p1 sig_p1],[0 500],'b--','LineWidth',2);xlabel('Process Var,\sigma^2_{q1}');ylabel('Frequency')
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));    
    saveas(gcf,'varHistQ1Nonlinear.jpg')
    figure;
    histogram(var_vect(3,:),'FaceColor','w','EdgeColor','k'),hold on;
    plot([sig_p2 sig_p2],[0 500],'b--','LineWidth',2);xlabel('Process Var,\sigma^2_{q2}');ylabel('Frequency')
    set(gca,'FontSize',16)
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'varHistQ2Nonlinear.jpg')
    % histograms of position estimate error (Mahalanobis)
    figure;
    histogram(distMah_irls,[0:2500:50000],'FaceColor','w','EdgeColor','k');hold on;
    xlabel('Mahalanobis Distance of State Est Error');ylabel('Frequency');
    plot([2*n 2*n],[0 1000],'--r','LineWidth',2)
    %axis([0 400 0 20]);
    set(gca,'FontSize',16);
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(gcf,'MahHistNonlinear.jpg');
    disp(['Mah Dist 2.5, 50, 97.5 Percentiles: ',num2str(round(prctile(distMah_irls,[2.5 50 97.5]),3))])
end