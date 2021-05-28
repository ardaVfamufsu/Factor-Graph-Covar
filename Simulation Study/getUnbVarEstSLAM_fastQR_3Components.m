%find unbiased estimate of process and measurement variances - use QR factorization result
function [sig_pHat1_std, sig_pHat2_std, sig_eHat_std]=getUnbVarEstSLAM_fastQR_3Components(R,rx1,rx2,rz,Ax1,Ax2,Az,flag_scale_cov_change)

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

% find traces of H1, H2, H3, H1^2, etc. without matrix inversion
M1=A1/R;
clear A1
M1_sq = M1*M1';
tH1=sum(diag(M1_sq));
tH11=sum(sum(M1_sq.*M1_sq));
M2=A2/R;
clear A2
M2_sq = M2 * M2';
tH2=sum(diag(M2_sq));
tH22=sum(sum(M2_sq.*M2_sq));
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

