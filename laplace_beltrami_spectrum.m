function [evecs evals error] = laplace_beltrami_spectrum(surf_filename, k)


try
    [V,F] = read_off(surf_filename);
    if size(V,2) > k
        % k = 50;
        V = normalize_by_SA(V, F);
        L = cotmatrix(V,F'); % some implementation may put L = 0.5*cotmatrix(V,F); this just scales ED by 0.5 and does not really matter.
        M = massmatrix(V,F','barycentric');
        [EV,ED] = eigs(L,M,k,'sm');
        ED = - ED; % minus sign is needed since cotmatrix() gives a negative definite matrix.

        evals = diag(ED);
        evecs = EV;
            error = 0;
    else 
        evals = [];
        evecs = [];
        error = 1;
    end
catch
    evals = [];
    evecs = [];
    error = 1;
end
end
