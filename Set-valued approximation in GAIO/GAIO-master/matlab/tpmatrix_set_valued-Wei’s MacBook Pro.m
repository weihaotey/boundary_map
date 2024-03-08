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
 % (C) 2021, Wei Hao Tey

 d = t.dim;
 b = t.boxes(depth); % get the geometry of the boxes
 N = size(b,2); S = whos('X'); l = floor(1e8/S.bytes/size(Xs,1));
 I = []; IJS = []; tic; length_each = zeros(floor(N/l)+2,1);
 if verbose, dispr('',0); end
 length_all = 0; mat_all = cell(1,floor(N/l)+1);
 for k = 0:floor(N/l) % split in chunks of size l
     t1 = t;
     K=k*l+1:min((k+1)*l,N);
     c = b(1:d,K); r = b(d+1:2*d,1); % center and radii of the boxes
     n = size(c,2); E = ones(n,1); % n is number of boxes in this loop
     P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); % sample points in all boxes
     P = f(P); % map the points
     ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
     eP = kron(Ei,Xs)*diag(epsilon) + kron(P,ones(size(Xs,1),1));

     %I = ones(length(eP),1);
     %for i = 1:length(eP)
     %   [minn,argmin]=min(vecnorm(eP(i,:)'-b(1:d,:),1));
     %   if minn < norm(b(d+1:2*d,1))
     %       I(i) = argmin;
     %   else
     %       I(i) = -1;
     %   end
     %   dispr(sprintf('%d \n',i,toc),1);
     %end

     I = t1.search(eP', depth); % get box numbers of image points
     pI = find(I>0); % I <= 0 iff f(P) not in state space
     J = kron(K',ones(size(X,1)*size(Xs,1),1)); % column numbers
     [I,J,S] = find(sparse(I(pI), J(pI), 1, N, N)); % transition matrix
     mat_all{k+1} = [I,J,S];
     %IJS = [IJS; I,J,S];
     length_all = length_all + length(I);
     length_each(k+2) = length_all;
     if (verbose)
        dispr(sprintf('%d of %d boxes, %.1f sec \n',min((k+1)*l,N),N,toc),1);
     end
 end
IJS = zeros(length_all,3);
for k = 0:floor(N/l)
    I = mat_all{k+1};
    IJS(length_each(k+1)+1:length_each(k+2),:) = I;
end
res = sparse(IJS(:,1), IJS(:,2), IJS(:,3), N, N); % transition matrix
cs = sum(res); % normalization generates stochastic matrix (transposed)
res = res*spdiags(1./cs', 0, N, N); % but probably changes dynamics
disp(toc)
if verbose, fprintf('\n'); end