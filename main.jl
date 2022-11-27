include("VTFM/src/VTFM.jl")
# include("VTFM/src/sumtest.jl")

using .VTFM
VTFM.greet()
ns = [2, 2, 2]
ds = [1, 1, 1]
s = ftest(1,1);
println("# The value of my var is: ", s)