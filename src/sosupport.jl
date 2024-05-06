#  ===============  sosupport.jl ======================================================================

#  Julia code to check whether a given SoS form f of degree 2k is uniquely Sum-of-Squares representable. 
#  Using semidefinite programming, two Gram matrix representations of f are computed: 
#  One lies in the relative interior of the Gram spectrahedron, the other one on the boundary. 
#  The two Gram matrices have the same rank if and only if Gram(f) is a singleton set, i.e. 
#  iff f is uniquely Sum-of-Squares representable. 
#  The code also computes whether the supporting face of f in Σ_2k is exposed. 

using DynamicPolynomials, SumOfSquares, LinearAlgebra;
using DataFrames;

import MultivariatePolynomials.degree
degree(f::Polynomial{T} where T) = max(degree.(monomials(f))...) # define degree of a polynomial

if !@isdefined default_solver
    using COSMO;
    default_solver = Symbol(COSMO);
end

mutable struct SOSDataContainer
    M_E;
    G;
    G_boundary;
    SOSDataContainer() = new();
end

mutable struct SOSData 
    f # sos polynomial
    solver;
    silent;
    k; # = degree(f) ÷ 2; 
    vars;
    n;
    digits;
    uniquely_sos_representable::Bool; # true iff f is uniquely sos representable 
    dual_nondegenerate::Bool; # true iff the supporting face of f in Σ is exposed. Equivalent: rk(G) == dim ker(M_E)
    dim_sosupp::Int; # rank of G
    sosupp;  # image of G
    # some relative interior point of the Gram spectrahedron Gram(f) 
    G; # Gram Matrix
    σ; # SOS decomposition of f corresponding to G
    λ; # eigenvalues of G
    monomials; # monomials used in relint sos representation
    # some boundary point of the Gram spectrahedron Gram(f) 
    G_boundary; # Matrix 
    σ_boundary; # SOS decomposition of f corresponding to G_boundary
    λ_boundary; # eigenvalues of G_boundary
    M_E; # moment matrix of relative interior point of dual cone C_f cut out by the hyperplane {E | E(f) = 0}.
    μ; # eigenvalues of M_E
    spectral_gap::SOSDataContainer; # spectral gaps of G, G_boundary and M_E. 
    rank::SOSDataContainer; # ranks of G, G_boundary and M_E. 
    corank::SOSDataContainer; # coranks of G, G_boundary and M_E. 
    solve_time; # solver time to compute G and M_E
    solve_time_boundary; # solver time to compute G_boundary
    SOSData(f, solver, silent) = new(f, Symbol(solver), silent, degree(f) ÷ 2, variables(f), length(variables(f)));
end

SOSData(f) = SOSData(f, default_solver, true);
spectrum(M) = eigen(Matrix(M)).values; 


calc_sos_attributes = function(f; digits=7, solver=default_solver, silent=true, feas_tolerance=1e-12, opt_tolerance=2000)
    s = SOSData(f, solver, silent);
    calc_sos_attributes!(s; digits=digits, feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance)
    return s;
end

calc_sos_attributes! = function(s::SOSData; digits=7, feas_tolerance=1e-12, opt_tolerance=2000)
    s.digits=digits;
    G, M_E, s.solve_time = relint_sos_representation(s.f; solver=s.solver, feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance);
    s.σ = SOSDecomposition(G, 10.0^(-digits-1));
    s.G = G.Q;
    s.M_E = M_E.Q;
    s.monomials = G.basis.monomials;
    s.λ = spectrum(s.G);
    s.μ = spectrum(s.M_E);
    G_boundary, s.solve_time_boundary = boundary_sos_representation(s.f, s.monomials; solver=s.solver, feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance);
    s.G_boundary = G_boundary.Q;
    s.σ_boundary = SOSDecomposition(G_boundary, 10.0^(-digits-1));
    s.λ_boundary = spectrum(s.G_boundary);
    calc_dimensions!(s, digits)
    calc_sosupp!(s);
end

calc_dimensions! = function(s::SOSData, digits)
    s.spectral_gap = SOSDataContainer();
    s.rank = SOSDataContainer();
    s.corank = SOSDataContainer();
    for (v, M) in zip([s.μ, s.λ, s.λ_boundary], [:M_E, :G, :G_boundary])
        v_rounded = round.(v, digits=digits);
        num_zeros = sum(v_rounded .== 0);
        num_nonzeros = sum(v_rounded .≠ 0);
        setfield!(s.rank, M, num_nonzeros);
        setfield!(s.corank, M, num_zeros);
        min_nontruncated = num_zeros < length(v) ? v[num_zeros + 1] : Inf; # min(v_rounded[v_rounded .≠ 0]...);
        max_truncated = num_zeros > 0 ? v[num_zeros] : -Inf; # max(abs.(v-v_rounded)...);
        setfield!(s.spectral_gap, M, (min_nontruncated, max_truncated));

    end
    s.dim_sosupp = s.rank.G;
    s.uniquely_sos_representable = (s.rank.G == s.rank.G_boundary);
    s.dual_nondegenerate = (s.rank.G == s.corank.M_E);
end

calc_sosupp! = function(s::SOSData)
    eigvecs = eigen(Matrix(s.G)).vectors[:, end-s.rank.G+1:end];
    s.sosupp = [v⋅s.monomials for v in eachcol(eigvecs)];
end

set_mosek_attributes = function(model; feas_tolerance=1e-16, opt_tolerance=1000)
    # confer: https://docs.mosek.com/10.0/capi/parameters.html#mosek.dparam.intpnt_co_tol_pfeas
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL", opt_tolerance); # default: 1000 # Optimality tolerance used by the interior-point optimizer for conic problems. If MOSEK cannot compute a solution that has the prescribed accuracy then it will check if the solution found satisfies the termination criteria with all tolerances multiplied by the value of this parameter. If yes, then the solution is also declared optimal.
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",   feas_tolerance); # Dual feasibility tolerance used by the interior-point optimizer for conic problems.
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS",   feas_tolerance); # Primal feasibility tolerance used by the interior-point optimizer for conic problems.
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_MU_RED",  feas_tolerance); # Relative complementarity gap tolerance used by the interior-point optimizer for conic problems.
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", feas_tolerance); # Relative gap termination tolerance used by the interior-point optimizer for conic problems.
end

relint_sos_representation = function(f₂; solver=default_solver, silent=true, feas_tolerance=1e-12, opt_tolerance=2000) 
    model = SOSModel(eval(solver).Optimizer);
    if silent; set_silent(model); end;
    if Symbol(solver)==:Mosek; set_mosek_attributes(model; feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance); end
    @constraint(model, c1, f₂ >= 0, maxdegree = degree(f₂))
    optimize!(model)
    if !silent
        @show termination_status(model);
        @show objective_value(model);
    end
    solve_time = MOI.get(model, MOI.SolveTimeSec())
    M_E = moment_matrix(c1);
    G = gram_matrix(c1);
    return G, M_E, solve_time;
end

boundary_sos_representation = function(f₂, monomials; solver=default_solver, silent=true, feas_tolerance=1e-12, opt_tolerance=2000) 
    model = Model(eval(solver).Optimizer);
    if silent; set_silent(model); end;
    if Symbol(solver)==:Mosek; set_mosek_attributes(model; feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance); end
    #monomials = vec([X[i]*Y[j] for i=1:n, j=1:n])
    N = length(monomials); 
    @variable(model, G[1:N, 1:N] in PSDCone())
    @constraint(model, f₂ == monomials'*G*monomials)
    A = Symmetric(randn(N, N));
    @objective(model, Max, tr(A*G))
    optimize!(model)
    if !silent
        @show termination_status(model)
        @show objective_value(model)
    end
    solve_time = MOI.get(model, MOI.SolveTimeSec())
    G_opt = value.(G);
    gram_G = GramMatrix(G_opt, MonomialBasis(monomials));
    return gram_G, solve_time;
end

is_uniquely_sos_representable = function(f₂, digits; solver=default_solver) 
    G1, M_E, _ = relint_sos_representation(f₂);
    G1 = Matrix(G1.Q); 
    M_E = Matrix(M_E.Q);
    G2, _ = boundary_sos_representation(f₂);
    λ1 = round.(eigen(G1).values, digits=digits);
    λ2 = round.(eigen(G2).values, digits=digits);
    μ = round.(eigen(M_E).values, digits=digits);
    rank_difference = sum(λ1 .≠ 0) - sum(λ2 .≠ 0);
    uniquely_sos_representable = (rank_difference == 0);
    dual_nondegenerate = (sum(λ1 .≠ 0) == sum(μ .== 0));
    return uniquely_sos_representable;
end

report = function(s::SOSData, m=nothing; transpose=false)
    # create a vector of field names
    fields = [:n, :k, :solver, :solve_time, :digits, :uniquely_sos_representable, :dual_nondegenerate, :dim_sosupp]
    special_fields = [:spectral_gap_G, :corank_dual, :spectral_gap_M_E, :rank_G_boundary, :spectral_gap_G_boundary];
    special_fields_values = [s.spectral_gap.G, s.corank.M_E, s.spectral_gap.M_E, s.rank.G_boundary, s.spectral_gap.G_boundary];
    explanations = ["No. of variables", "degree", "SDP solver", "SDP solve time for (G, M_E)", "Threshold for truncating eigenvalues (in digits)", "", "Whether f lies on exposed face of Σ_2k", "", "(smallest nontruncated, largest truncated) eigval", "corank of M_E", "(smallest nontruncated, largest truncated) eigval", "", ""]

    # create data frame with existing parameters and fields of s
    if m === nothing
        df = DataFrame(
            Parameter = vcat(fields, special_fields),
            Value = vcat([getfield(s, f) for f in fields], special_fields_values),
            Explanation = explanations
        )
    else
        df = DataFrame(
            Parameter = vcat([:m], fields, special_fields),
            Value = vcat([m], [getfield(s, f) for f in fields], special_fields_values),
            Explanation = vcat(["No. of addends"], explanations)
        )
    end
    if transpose
        return report_transpose(df)
    end
    return df;
    #(:f, :solver, :silent, :digits, :σ, :σ_boundary, :uniquely_sos_representable, :dual_nondegenerate, :dim_sosupp, :sosupp, :G, :monomials, :G_boundary, :λ, :λ_boundary, :M_E, :μ, :spectral_gap, :rank, :corank)
end


# returns a DataFrame row containing just the data of the report
report_transpose = function(df::DataFrame)
    df_result = empty_sosreport(write_m=:m in df.Parameter);
    push!(df_result, df.Value);
    return df_result;
    # df_tp = permutedims(df);
    # # set the first row as the header of the transposed table
    # rename!(df_tp, pairs(df_tp[1,:]) |> Dict);
    # # select only the row with data from the transposed table
    # return df_tp[2:2, :];
end



# reportall = function(S::Vector{SOSData})
#     fields = [:n, :k, :solver, :digits, :uniquely_sos_representable, :dual_nondegenerate, :dim_sosupp];
#     special_fields = [:spectral_gap_G, :corank_dual, :spectral_gap_M_E, :rank_G_boundary, :spectral_gap_G_boundary];    
#     df = DataFrame([field => [] for field in vcat(fields, special_fields)])
#     for s in S
#         df_s = report(s);
#         df_new_row = report_transpose(df_s)
#         append!(df, [df_new_row]);
#     end
#     return df;
# end


empty_sosreport = function(; write_m=false)
    fields = [:n, :k, :solver, :solve_time, :digits, :uniquely_sos_representable, :dual_nondegenerate, :dim_sosupp];
    types = [Int, Int, Symbol, Float64, Int, Bool, Bool, Int];
    special_fields = [:spectral_gap_G, :corank_dual, :spectral_gap_M_E, :rank_G_boundary, :spectral_gap_G_boundary];
    special_types = [Tuple{Float64, Float64}, Int, Tuple{Float64, Float64}, Int, Tuple{Float64, Float64}];
    if write_m
        df = DataFrame([field => type[] for (field, type) in zip(vcat([:m], fields, special_fields), vcat([Int], types, special_types))]);
    else
        df = DataFrame([field => type[] for (field, type) in zip(vcat(fields, special_fields), vcat(types, special_types))]);
    end
    return df;
end
