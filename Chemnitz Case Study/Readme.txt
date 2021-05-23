File dependencies:
main_ChemnitzMestimateSagnacIRLS_3component
	Data_Chemnitz.csv

Code usage

flag_unbiased=1;
flag_scale_cov_change=1;
[x_sol, rhs_x, rhs_z, Ax, Az, A,rhs,sig_pHat1,sig_pHat2, sig_eHat]=main_ChemnitzMestimateSagnacIRLS_3component(flag_unbiased,flag_scale_cov_change);
