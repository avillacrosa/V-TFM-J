module VTFM
    include("meshgen.jl")
    include("runFEM.jl")
    include("internalF.jl")
    include("stress.jl")
    include("deformF.jl")
    include("getdNdx.jl")
    include("getB.jl")
    include("constK.jl")
    include("ctensor.jl")
    include("BCtoNodal.jl")
    include("assembleK.jl")
    include("constKSparse.jl")

    export constKSparse
    export assembleK
    export ctensor
    export meshgen
    export runFEM
    export internalF
    export stress
    export deformF
    export getdNdx
    export getB
    export constK
    export BCtoNodal
end
