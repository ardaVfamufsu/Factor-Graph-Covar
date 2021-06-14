n_tst = 2000;
clear samp
samp = zeros(3,n_tst);
for ii=1:n_tst
    tst = H *randn(120,1);
    for jj=1:3
        samp(jj,ii)=sum(tst(idx{jj}).^2);
    end
end

for ii=1:3
    figure(ii)
    hold off
    handle = histogram(samp(ii,:),'normalization','pdf');
    my_x = handle.BinLimits(1):.5:handle.BinLimits(2);
    my_y = chi2pdf(my_x,sum(T(ii,:)));
    hold on;
    plot(my_x,my_y)
end

% create chisquare probability plot of the data
% if the points plot near the straight line then the data approximately
% follows the chisquare distribution
% chisquare distribution is just a gamma with alpha parameter nu/2 (nu = dof) and beta parameter = 2
df=length(idx{1}); % all 3 residual vectors have same lengths
pd = makedist('Gamma','a',df/2,'b',2);
figure
subplot(1,3,1),probplot(pd,samp(1,:));title('process resid 1');
subplot(1,3,2),probplot(pd,samp(2,:));title('process resid 2');
subplot(1,3,3),probplot(pd,samp(3,:));title('measure resid 1');
% Test the null hypothesis that the data comes from a population with a chisquare distribution. 
% Use fitdist to create a probability distribution object with parameters estimated from the data.
% if h =1 (or p is small, e.g., <0.05) then the null hypothesis that chi-sq distribution fits the data well is rejected
pd1 = fitdist(samp(1,:)','Gamma');
[h1 p1] = chi2gof(samp(1,:),'CDF',pd1)
pd2 = fitdist(samp(2,:)','Gamma');
[h2 p2] = chi2gof(samp(2,:),'CDF',pd2)
pd3 = fitdist(samp(3,:)','Gamma');
[h3 p3] = chi2gof(samp(3,:),'CDF',pd3)

%% Conclusion:  pdfs and histograms do not match
%% 3 is always under.  1&2 are always over
%% The sum(T(ii,:)) is much closer than if I add up the 
%% pdfs externally though!
