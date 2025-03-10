function res = tpmatrix_set_valued(t, f, epsilon, X, Xs, depth, verbose,multiplicative)
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
 % multiplicative: '0' or '1' ('1' use multiplicative noise)
 % (C) 2013, djs GbR
 % (C) 2015, Hans Peschke
 % (C) 2021, Wei Hao Tey

 d = t.dim;
 fast = 0; % whether to have fast computation
 b = t.boxes(depth); % get the geometry of the boxes
 N = size(b,2); S = whos('X'); l = floor(5e8/S.bytes/size(Xs,1));
 I = []; IJS = []; tic;
 if verbose, dispr('',0); end
 if fast, length_all = 0; mat_all = cell(1,floor(N/l)+1); end
 for k = 0:floor(N/l) % split in chunks of size l
     K=k*l+1:min((k+1)*l,N);
     c = b(1:d,K); r = b(d+1:2*d,1); % center and radii of the boxes
     n = size(c,2); E = ones(n,1); % n is number of boxes in this loop
     P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); % sample points in all boxes
     oldP = P; % make a copy of P for changing bound
     P = f(P); % map the points
%      P = f(P); % map the points again
     ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
     if multiplicative, Ei = vecnorm(P-oldP,2,2);end % multiplicative noise
     eP = kron(Ei,Xs)*diag(epsilon) + kron(P,ones(size(Xs,1),1));
     I = t.search(eP', depth); % get box numbers of image points
     pI = find(I>0); % I <= 0 iff f(P) not in state space
     J = kron(K',ones(size(X,1)*size(Xs,1),1)); % column numbers
     [I,J,S] = find(sparse(I(pI), J(pI), 1, N, N)); % transition matrix
     clear eP
     if fast
        mat_all{k+1} = [I,J,S]; % record all the splitted martrix value
        length_all = length_all + length(I); % record total length of matrix value
        if length_all > 1e10
            disp(length_all)
            kk = zeros(length_all,3);
            error('memory insufficient')
        end
     else
         IJS = [IJS; I,J,S];
     end
     if (verbose)
        dispr(sprintf('%d of %d boxes, %.1f sec \n',min((k+1)*l,N),N,toc),1);
     end
 end
 if fast
     if verbose, fprintf('%d x %d matrix \n',length_all,3); end
    % initialise vector for position index IJ and transition probability S
    IJ = zeros(length_all,2,'int32'); count = 0; S = zeros(length_all,1);
    % put all the splitted part to form matrix
    for k = 0:floor(N/l)
        I = mat_all{k+1};
        IJ(count+1:count+length(I),:) = I(:,1:2);
        S(count+1:count+length(I)) = I(:,3);
        count = count + length(I);
    end
    clear I
    res = sparse(IJ(:,1), IJ(:,2), S, N, N); % transition matrix
 else
    res = sparse(IJS(:,1), IJS(:,2), IJS(:,3), N, N); % transition matrix
 end
% cs = sum(res); % normalization generates stochastic matrix (transposed)
res = res*spdiags(1./sum(res)', 0, N, N); % but probably changes dynamics

if verbose, fprintf('%f seconds\n',toc); end