Pt = tpmatrix(t, f, eps, .) % compute the transpose of the transition
% probability matrix P with the Ulam-method
[Vl, ~] = eigs(Pt,1,'lr'); % compute the eigenvector of the largest
% real eigenvalue of the transpose of P
M = (Vl > c); % choose the boxes with strict positive
% probability (c > 0)