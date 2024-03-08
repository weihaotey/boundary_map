function res = tpmatrix_set_valued_inverse(t, f, eps, X, Xs, depth, verbose)
 % TPMATRIX transition probability matrix with epsilon inflation
 %
 % TPMATRIX(t, f, p, eps, X, Xs, d, v) computes the matrix of
 % transition probabilities between the boxes in the tree t
 % t GAIO-tree containing the box covering
 % f single-valued transformation
 % eps epsilon for inflation the image
 % X m x d-matrix of sample points
 % Xs matrix of sample points for the inflated image
 % d depth of the tree on which the matrix is computed
 % v verbose flag: '0' or '1' ('1' prints progress)
 % (C) 2013, djs GbR
 % (C) 2015, Hans Peschke

 d = t.dim;
 b = t.boxes(depth); % get the geometry of the boxes
 N = size(b,2); S = whos('X'); l = floor(5e7/S.bytes);
 I = []; IJS = []; tic;
 if verbose, dispr('',0); end
 for k = 0:floor(N/l) % split in chunks of size l
 K=k*l+1:min((k+1)*l,N);
 c = b(1:d,K); r = b(d+1:2*d,1); % center and radii of the boxes
 n = size(c,2); E = ones(n,1); % n is number of boxes in this loop
 P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); % sample points in all boxes
 ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
 eP = kron(Ei,Xs)*diag(eps) + kron(P,ones(size(Xs,1),1));
 eP = f(eP); % map the points
 I = t.search(eP', depth); % get box numbers of image points
 pI = find(I>0); % I <= 0 iff f(P) not in state space
 J = kron(K',ones(size(X,1)*size(Xs,1),1)); % column numbers
 [I,J,S] = find(sparse(I(pI), J(pI), 1, N, N)); % transition matrix
 IJS = [IJS; I,J,S];
 if (verbose)
 dispr(sprintf('%d of %d boxes, %.1f sec',min((k+1)*l,N),N,toc),1);
 end
end
res = sparse(IJS(:,1), IJS(:,2), IJS(:,3), N, N); % transition matrix
cs = sum(res); % normalization generates stochastic matrix (transposed)
res = res*spdiags(1./cs', 0, N, N); % but probably changes dynamics
if verbose, fprintf('\n'); end