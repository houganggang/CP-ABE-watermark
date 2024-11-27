function [Vnew] = normalize_by_SA(V, F)

    SA = trimeshSurfaceArea(V', F');
    SqrtSA = sqrt(SA);
    Vnew = V'./SqrtSA;
end
    %SA = trimeshSurfaceArea(Vnew, F');
    



