function res = tpmatrix_set_valued(t, f, epsilon, X, Xs, depth, verbose)
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
 P = f(P); % map the points
 ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
 eP = kron(Ei,Xs)*diag(epsilon) + kron(P,ones(size(Xs,1),1));
 I = t.search(eP', depth); % get box numbers of image points
 J = kron(K',ones(size(X,1)*size(Xs,1),1)); % column numbers
 J = J(I>0); I = I(I>0); 
 pt = unique([I,J],'rows'); I = pt(:,1); J = pt(:,2);
 [I,J,S] = find(sparse(I, J, 1, N, N)); % transition matrix
 IJS = [IJS; I,J,S];
 if (verbose)
 dispr(sprintf('%d of %d boxes, %.1f sec \n',min((k+1)*l,N),N,toc),1);
 end
end
res = sparse(IJS(:,1), IJS(:,2), IJS(:,3), N, N); % transition matrix
if verbose, fprintf('\n'); end
%{
 d = t.dim;
 b = t.boxes(depth); % get the geometry of the boxes
 N = size(b,2); S = whos('X'); l = floor(5e7/S.bytes);
 tic; II = []; JJ = [];
 for k = 1:N % split in chunks of size l
 c = b(1:d,k); r = b(d+1:2*d,1); % center and radii of the boxes
 n = 1; E = ones(n,1); % n is number of boxes in this loop
 P = f(c'); % map the points
 ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
 eP = Xs*diag(epsilon) + kron(P,ones(size(Xs,1),1));
 I = unique(t.search(eP', depth)); % get box numbers of image points
 I = I(I>0); JJ = [JJ;I]; II = [II;k*ones(length(I),1)];
 end
res = sparse(II,JJ,1,N,N); %Binary matrix
toc;
 %}