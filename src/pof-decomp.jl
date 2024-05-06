#  ===============  pof-decomp.jl ======================================================================

#  Contains methods to decompose some third-order powers-of-forms decompositions. 
#  The methods work e.g. when the second order power sum is uniquely sos representable. 

include("sosupport.jl");


mutable struct PofDecomposition 
    f_2; # sos polynomial of degree 2*k
    f_3; # positive weighted sum of cubes of k-forms 
    k; # degree of addends
    vars;
    n; # == length(vars)
    q; # addends of pof decomposition
    λ; # weights of pof decomposition
    sosdata; # SOSData object associated with the Sum-Of-Squares decomposition of f_2 from sosupport.jl 
    digits;
    u; # generating system of some superspace U of <q_1,…,q_m>
    N; # == length(u)
    Y; # new variables Y_1,…,Y_N
    M_φ1;
    M_φ2; # representing matrix of φ: Y_i ↦ u_i on graded components 1, 2, 3. 
    M_φ3; 
    ℓ;   # == φ^{-1}(q)
    g_2; # == φ^{-1}(f_2)
    g_3; # == φ^{-1}(f_3)
    PofDecomposition(f_2, f_3) = new(f_2, f_3, degree(f_3) ÷ 3, variables(f_2 + f_3), length(variables(f_2 + f_3)));
end


raw"""
    Semidefinite algorithm for powers-of-forms decomposition. 
    Inputs: Polynomials 
    + ``f_2``: of degree ``2k``.
    + ``f_3``: of degree ``3k``. 
    for some ``k \in \mathbb N``. 

    Assumption: There exist k-forms ``q_1,\ldots,q_m``  and positive reals $\lambda_1,\ldots,\lambda_m$ such that 
    + ``f_2 = \sum_{i=1}^m \lambda_i q_i^2``.
    + ``f_3 = \sum_{i=1}^m \lambda_i q_i^3``. 

    Outputs: 
    + ``(\lambda_1,q_1),\ldots,(\lambda_m,q_m)``, in any order. 
"""
pof_decompose = function(f_2, f_3; digits=7, u=nothing)
    pof = PofDecomposition(f_2, f_3);
    if degree(f_2) ÷ 2 ≠ pof.k 
        error("Degrees of f_2 and f_3 are not compatible.");
        return;
    end
    
    if u === nothing
        # if no space was given, calculate the sum of squares support. 
        # parameters are set to high accuracy. 
        pof.sosdata = calc_sos_attributes(f_2, digits=digits, feas_tolerance=1e-16, opt_tolerance=8000);
        pof.u = pof.sosdata.sosupp;
    else 
        pof.u = u;
        pof.N = length(u);
    end
    pof.N = length(pof.u);
    μ_min = set_φ_matrices!(pof);
    is_independent = round.(μ_min; digits=digits) > 0;
    if !is_independent 
        error("Basis of the space U has algebraic dependencies of degree 3!\n Minimum eigenvalue, treated as zero: $(μ_min)")
        return;
    end
    @polyvar Y[1:pof.N];
    pof.Y = Y;
    compute_φ_preimages!(pof, f_2, f_3);
    pof.ℓ, pof.λ = positive_weighted_sylvester(pof.g_2, pof.g_3);
    pof.q = [l(pof.Y=>pof.u) for l in pof.ℓ]; # substitute back via φ
    return pof;
end

is_cubically_independent = function(u; digits=7)
    N = length(u);
    k = degree.(u[1].x[1])
    vars = variables(u);
    monoms = monomials(vars, 3k);
    UUU = [u[i]*u[j]*u[k] for i=1:N, j=1:N, k=1:N if i≤j≤k];
    M = hcat([coefficients(p, monoms) for p in UUU]...);
    μ_min = min(svd(M).S...); # minimum singular value
    is_independent = round.(μ_min; digits=digits) > 0;
    return is_independent, μ_min;
end
    
set_φ_matrices! = function(pof::PofDecomposition)
    u = pof.u;
    k = pof.k;
    monoms_1 = monomials(pof.vars, k);
    monoms_2 = monomials(pof.vars, 2k);
    monoms_3 = monomials(pof.vars, 3k);
    pof.M_φ1 = hcat([coefficients(u[i], monoms_1) for i=1:pof.N]...);
    pof.M_φ2 = hcat([coefficients(u[i]*u[j], monoms_2) for i=1:pof.N for j=1:pof.N if i<=j]...);
    pof.M_φ3 = hcat([coefficients(u[i]*u[j]*u[k], monoms_3) for i=1:pof.N for j=1:pof.N for k=1:pof.N if i<=j<=k]...);
    μ_min = min(svd(pof.M_φ3).S...); # check minimum singular value to see if u_iu_ju_k are linearly dependent
    return μ_min;
end

compute_φ_preimages! = function(pof, f_2, f_3)
    monoms_2 = monomials(pof.vars, 2pof.k)
    monoms_3 = monomials(pof.vars, 3pof.k);
    f_2_vec = coefficients(f_2, monoms_2);
    f_3_vec = coefficients(f_3, monoms_3);
    pof.g_2 = pof.M_φ2\f_2_vec ⋅ monomials(pof.Y, 2);
    pof.g_3 = pof.M_φ3\f_3_vec ⋅ monomials(pof.Y, 3);
    pof.g_2 = pof.g_2(pof.Y=>reverse(pof.Y)); # the monomials function by default sorts Y as Y[N],...,Y[1]
    pof.g_3 = pof.g_3(pof.Y=>reverse(pof.Y)); # here, we correct for this
end

positive_weighted_sylvester = function(g_2, g_3)
    Y = variables(g_2);
    v = randn(length(Y));
    g_v = (1.0/3)*differentiate(g_3, Y)⋅v;
    u, M_2, M_v = sylvester_matrix_representations(g_2, g_v);
    m = length(u); # rank of g_2. Note: might be smaller than the number of variables Y!
    eig_dec = eigen(M_v, M_2); # generalized eigendecomposition
    μ = eig_dec.values;
    b = [vec for vec in eachcol(Symmetric(M_2)*eig_dec.vectors)]; # dual of generalized Eigenvectors == rk1 terms  (up to multiples)
    u_v = [u[i](Y=>v) for i=1:m]; # evaluate u in v. 
    a = [(μ[j]/(b[j]⋅u_v))*b[j] for j=1:m]; # correct the multiples 
    ℓ = [a[i]⋅u for i=1:m]; # convert back to linear forms
    λ = sylvester_find_weights(g_2, ℓ); # find the missing weights
    return ℓ, λ;
end

"""
    sylvester_matrix_representations

    input: 
    +  g_2, g_v: two quadratic forms which have the same image.
    output: 
    1. u: a basis of the image space of g_2, represented by polynomials.  
    2. M_2: some matrix such that g_2 = u^T * M_2 * u. 
    3. M_v: some matrix such that g_v = u^T * M_v * u.
"""
sylvester_matrix_representations = function (g_2, g_v; digits=7)
    Y = variables(g_2);
    N = length(Y);
    qf_to_symmatrix(qf) = [(i==j ? 1 : 1/2) * float(DynamicPolynomials.coefficient(qf, Y[i]*Y[j])) for i=1:N, j=1:N];

    G_2 = qf_to_symmatrix(g_2);
    G_v = qf_to_symmatrix(g_v);

    eig_dec = eigen(G_2);
    eig_vals = eig_dec.values;
    
    nonzero_index = min([i for i=1:N if round(eig_vals[i]; digits=7) > 0]...);
    if nonzero_index == 1 # if g_2 has full rank, nothing needs to be done. 
        return Y, G_2, G_v;
    end
    B = eig_dec.vectors[:, nonzero_index:end];
    M_2 = B'*G_2*B;
    M_v = B'*G_v*B;
    u = [vec⋅Y for vec in eachcol(B)];
    return u, M_2, M_v; 
end


sylvester_find_weights = function (g_2, ℓ)
    Y = variables(ℓ); 
    m = length(ℓ);
    monoms_2 = monomials(Y, 2);
    L = hcat([coefficients(ℓ[i]^2, monoms_2) for i=1:m]...);
    g_2_vec = coefficients(g_2, monoms_2);
    μ = L\g_2_vec;
    return μ;
end

sylvester_produce_g2 = function(g_3)
    Y = variables(g_3); 
    m = length(Y);
    qf_to_symmatrix(qf) = [(i==j ? 1 : 1/2) * float(DynamicPolynomials.coefficient(qf, Y[i]*Y[j])) for i=1:m, j=1:m];
    M = [];
    for i=1:m
        g_i = (1.0/3).*differentiate(g_3, Y[i]);
        push!(M, qf_to_symmatrix(g_i));
    end
    model = Model(eval(solver).Optimizer);
    @variable(model, a[1:m]);
    @constraint(model, sum(a[i]*M[i] for i=1:m) in PSDCone());
    optimize!(model);
    return Y'*sum(value(a[i])*M[i] for i=1:m)*Y
end


polycompare = function (p_1, p_2)
    if isless(leadingmonomial(p_1), leadingmonomial(p_2))
        return true;
    elseif isequal(leadingmonomial(p_1), leadingmonomial(p_2))
        m = leadingmonomial(p_1);
        c_1 = DynamicPolynomials.coefficient(p_1, m);
        c_2 = DynamicPolynomials.coefficient(p_2, m);
        if c_1 < c_2
            return true;
        elseif c_1 == c_2
            return polycompare(p_1 - c_1*m, p_2 - c_1*m);
        end
    end 
    return false;
end