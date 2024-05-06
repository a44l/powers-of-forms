# this file verifies that the family q = [(X[i] + X[j] + X[k])^2 for i=1:n, j=1:n, k=1:n if (i+j+k)%n == 0 && i≤j≤k] is cubically independent
# we checked cubic independence in the range n<=15. 
# notice that this also shows cubic independence of  q_2 = [(X[i] + X[j] + X[k])(Y[i] + Y[j] + Y[k]) for i=1:n, j=1:n, k=1:n if (i+j+k)%n == 0 && i≤j≤k],
# since the substitution Y=>X maps q_2 to q. 

using DynamicPolynomials, LinearAlgebra;

is_cubically_independent = function(u; digits=7)
    N = length(u);
    vars = variables(u);
    deg = degree.(u[1].x[1])
    monoms = monomials(vars, 3deg);
    UUU = [u[i]*u[j]*u[k] for i=1:N, j=1:N, k=1:N if i≤j≤k];
    M = hcat([coefficients(p, monoms) for p in UUU]...);
    μ_min = min(svd(M).S...); # minimum singular value
    is_independent = round.(μ_min; digits=digits) > 0;
    return is_independent, μ_min;
end

n = 9;
@polyvar X[1:n];
q = [(X[i] + X[j] + X[k])^2 for i=1:n, j=1:n, k=1:n if (i+j+k)%n == 0 && i≤j≤k] 
m = length(q);
is_cubically_independent(q)
    