function v = invariant_box(t, f, eps, X, Xs, depth)
 d = t.dim;
 b = t.boxes(depth); % get the geometry of the boxes
 N = size(b,2);
 v = zeros(N,1);
 for i=1:N
   c = b(1:d,i); r = b(d+1:2*d,1); % center and radii of the boxes
   n = size(c,2); E = ones(n,1); % n is number of boxes in this loop
   P = X*diag(r) + kron(c',ones(size(X,1),1)); % sample points in all boxes
   P = f(P); % map the points
   ni = size(P,1); Ei = ones(ni,1); % inflate the image points by epsilon
   eP = kron(Ei,Xs)*diag(eps) + kron(P,ones(size(Xs,1),1));
   I = t.search(eP', depth); % get box numbers of image points
   p = find(I==i);
   v(i) = length(p)/(size(X,1)*size(Xs,1));
 end
end