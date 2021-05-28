% find coordinates of 95% confidence interval from covariance Sig_x and
% estimates mu_x.
% measurement error ellipse coordinates: X,Y
function [X,Y,eigvec,eigval]=findConfEllipse(Sig_x,mu_x)
% unit circle
rc=1;xc=0;yc=0;th = 0:pi/50:2*pi;
xunit = rc * cos(th) + xc; yunit = rc * sin(th) + yc;
% confidence ellipse
[eigvec,eigval] = eig(Sig_x);
XY = [xunit' yunit']*sqrt(eigval)*eigvec';
r=2; %number of dimensions (degrees of freedom)
conf=0.95; k = sqrt(chi2inv(conf,r));
X = k*XY(:,1)+mu_x(1);
Y = k*XY(:,2)+mu_x(2);

end
