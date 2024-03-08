function [v, lambda, it, err] = power_method(A, v0, tol, maxit)
    it = 0; v = v0 / norm(v0);

    while(it == 0 || (max(abs(v - v_old)) > tol*lambda && it < maxit))
        v_old = v; it = it+1;
        v = (A*v); % eigenvector
        lambda = norm(v); % eigenvalue
        v = v / norm(v); % normalize
    end
    err = max(abs(v-v_old));
    v = v / max(v);
return