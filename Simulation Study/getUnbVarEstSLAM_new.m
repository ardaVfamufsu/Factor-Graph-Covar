function cov_est = getUnbVarEstSLAM_new(A,W,R,idx)

%% This function takes in an A, W, and R matrix (all unweighted) and
%% computes what the estimated covariance should be for different sets
%% of values.  The sets are specified by idx, where idx is a cell array
%% of indicies.  For example, idx{1}=1:10, idx{2}=11:20, and idx{3}=21:40.
%% Note that while I don't test it, idx should (a) cover the whole set
%% of rows of A (and W and R) and (b) each idx should be mutually exclusive
%% of the otherse.  In other words, if idx{1}=1:10, nothing between 1 and 10
%% should be in any of the other indices
%%
%% Args:
%%  A -- the (negative?) derivative of R w.r.t. the state values.  Shoudl be size mXn
%%  W -- the weighting matrix, size mXm (sparse?)
%%  R -- the residuals size mX1
%%  idx -- the division of the m indicies into sets
%% Returns:
%%  cov_est -- a vector of the covariance estimates.  indicies correspond
%%      with what is passed in by idx

%% Maybe I can test idx meets the conditions decribed above after all?
for ii = 1:length(idx)-1
    for jj = ii+1:length(idx)
        assert(isempty(intersection(idx{ii},idx{jj})))
    end
end

m = size(A,1);
full_idx = idx{1};
for ii=2:length(idx)
    full_idx=union(full_idx,idx{i});
end
assert(isempty(setdiff(full_idx,1:m)))
assert(isempty(setdiff(1:m,full_idx)))

%%Let's do some actual work here.
% First, form H.  No prettiness here, just brute force (i.e. may not work
% well with large graphs)
H = eye(m)-(A/(A.'*W*A))*A.'*W;

%Form the "T" matrix
T = zeros(length(idx));
Rsq = zeros(length(idx),1);
for ii = 1:length(idx)
    for jj=1:length(idx)
        small_H = H(idx{ii},idx{jj});
        T(ii,jj) = norm(small_H,'fro')^2;
    end
    Rsq = norm(R(idx{ii}),2)^2;
end

cov_est = T\Rsq;

    