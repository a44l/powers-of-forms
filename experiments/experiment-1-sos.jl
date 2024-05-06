using Pkg;
Pkg.activate("SoSupport");
using MosekTools; default_solver=Mosek; # change to your favourite solver. MOSEK needs a licence. Uncomment to use the COSMO solver as a default.  
include("../src/sosupport.jl");
include("../src/random-quadratics.jl");
using LinearAlgebra # e.g. for rank 
using DataFrames; using CSV; # to protocol results
using JLD2; using FileIO;

N = 10;
k = 2;
ranktol = 5; # measured in digits
feas_tolerance = 1e-8;
opt_tolerance = 1000;
@polyvar X[1:N] Y[1:N];

# -------------- paths to write results of computations
filepath_csv = "data/experiment-1-explicit-family.csv"; 
filepath_jld = "data/experiment-1-explicit-family.jld2"; 

# -------------- check if previously computed data exists 
restart = false; # set this to true in order to overwrite results of previous computations and start from scratch
df = empty_sosreport(write_m=true);
if isfile(filepath_csv)
    df = DataFrame(CSV.read(filepath_csv, DataFrame))
    if restart || (size(df, 1) == 0)
        df = empty_sosreport(write_m=true);
        CSV.write(filepath_csv, df);
        jldsave(filepath_jld);
    end
else 
    CSV.write(filepath_csv, df);
    jldsave(filepath_jld);
end
is_computed(n::Int) = 2n∈df.n


# -------------- 
for n = 2:N;
    if is_computed(n)
        continue;
    end
    q = [(X[i] + X[j] + X[k])*(Y[i] + Y[j] + Y[k])  for i=1:n, j=1:n, k=1:n if (i+j+k)%n == 0 && i≤j≤k] 
    m = length(q);
    
    # uncomment to check unique SOS representability in a neighbourhood
    # ε = 0.01;
    # q_noisy = q .+ ε.* [gaussian_splitvar_quadratic(X[1:n], Y[1:n]) for i=1:m]
    # f₂ = sum(q_noisy[i]^2 for i=1:m)
    
    f₂ = sum(q[i]^2 for i=1:m);
    
    s = calc_sos_attributes(f₂; digits=ranktol, feas_tolerance=feas_tolerance, opt_tolerance=opt_tolerance);
    df_new_row = report(s, m; transpose=true);
    append!(df, df_new_row; promote=true); # write to CSV
    # write to jld2
    jld_file = jldopen(filepath_jld, "a")
    write(jld_file, "s[$(n)]", s);
    close(jld_file);
    CSV.write(filepath_csv, df_new_row; append=true);
end




