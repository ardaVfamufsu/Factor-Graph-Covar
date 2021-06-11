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

%% Conclusion:  pdfs and histograms do not match
%% 3 is always under.  1&2 are always over
%% The sum(T(ii,:)) is much closer than if I add up the 
%% pdfs externally though!