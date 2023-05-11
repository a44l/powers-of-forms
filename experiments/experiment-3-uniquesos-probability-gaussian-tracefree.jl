using Pkg;
Pkg.activate("SoSupport");
using MosekTools; default_solver=Mosek;
include("../src/sosupport.jl");
include("../src/random-quadratics.jl");
using LinearAlgebra 
using DataFrames; using CSV; 
using JLD2; using FileIO;
using StatsBase;

N_min = 2;
N = 14;
ranktol = 5; # measured in digits
feas_tolerance = 1e-10;
opt_tolerance = 1000;
@polyvar X[1:N];

# -------------- paths to write results of computations
filepath = "data/experiment-3-gaussian-tracefree-study/"; 
seriesname = "m=n/";
M(n) = n;

# for each n=1:N, we generate R samples of a sum of M(n) Gaussian trace-free quadratics.
# we check for each whether their sum of squares is uniquely representable.  
R = 20;
for n=N_min:N 
    filepath_n_csv = string(filepath, seriesname, "samples-$n.csv");
    filepath_n_jld = string(filepath, seriesname, "samples-$n.jld2");
    df = empty_sosreport(write_m=true);
    if isfile(filepath_n_csv)
        continue;
    else
        mkpath(dirname(filepath_n_csv));
        CSV.write(filepath_n_csv, df);
        jldsave(filepath_n_jld);
    end
    m = M(n);
    for i=1:R
        q = [gaussian_trace0_quadratic(X[1:n]) for i=1:m]; 
        f₂ = sum(q[i]^2 for i=1:m);
        s = calc_sos_attributes(f₂; digits=ranktol, feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance);
        df_new_row = report(s, m; transpose=true);
        CSV.write(filepath_n_csv, df_new_row; append=true);
        append!(df, df_new_row; promote=true);
        jld_file = jldopen(filepath_n_jld, "a+");
        write(jld_file, "s[$i]", s);
        close(jld_file);
    end

end

# after computing a sample set for each n=1:N, we can compute the probability that f₂ is uniquely representable by averaging over the sample set. 
# we also compute the probability p_pointed for the stronger property that the sum of squares cone  
p_pointed = -ones(N);
p_unique_sos = -ones(N);

# outputs the generic Waring rank for squares of quadratic forms, as conjectured by Ottaviani.  
generic_rank = function(n)
    M_2 = binomial(n+1, 2);   # dimension of space of quadratics
    M_4 = binomial(n+3, 4); # dimension of space of quartics
    for m=M_2:-1:1 
        if m*M_2 - binomial(m, 2) < M_4
            return m + 1;
        end
    end
end

for n=N_min:N
    if M(n) > generic_rank(n)
        continue; # ignore cases where the forms q[i] are linearly dependent
    end
    df_n = DataFrame(CSV.read(string(filepath, seriesname, "samples-$n.csv"), DataFrame));
    p_pointed[n] = mean(df_n.m .== df_n.corank_dual);
    p_unique_sos[n] = mean(df_n.m .== df_n.dim_sosupp);
end
filepath_probabilities = string(filepath, seriesname, "probabilities.csv");
df_prob = DataFrame(n=1:N, p_pointed=p_pointed, p_unique_sos=p_unique_sos);
CSV.write(filepath_probabilities, df_prob);
