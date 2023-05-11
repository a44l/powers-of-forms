using Random;

gaussian_quadratic = function(variables)
    n = length(variables);
    A = randn(n,n);
    return variables'*A*variables
end;

gaussian_splitvar_quadratic = function(X, Y)
    n_1 = length(X);
    n_2 = length(Y);
    A = [zeros(n_1,n_1) randn(n_1, n_2); 
         randn(n_2,n_1) zeros(n_2, n_2)];
    return [X;Y]'*A*[X;Y]
end;

gaussian_trace0_quadratic = function(variables)
    n = length(variables);
    A = randn(n,n);
    A = A-1/n*tr(A)*I;
    return variables'*A*variables
end;

gauss_gramian_trace0_quadratic = function(variables)
    n = length(variables);
    A = 1/n .* randn(n,n);
    B = A'*A
    B = B-1/n*tr(B)*I;
    return variables'*B*variables
end;