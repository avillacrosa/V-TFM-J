include("VTFM/src/VTFM.jl")
using .VTFM
println("Running!")
@time begin
n = 5;
ns = [n, n, n]
ds = [1/(n-1), 1/(n-1), 1/(n-1)]
uBC = [3 0 3 0; 3 0 1 0; 3 0 2 0; 3 1 3 0.1]
E  = 100;
nu = 0.3;

runFEM(ns, ds, E, nu, uBC)

end
