module Ising

import LoadLeveller
import HDF5

mutable struct MC <: LoadLeveller.AbstractMC
    T::Float64

    spins::Array{Int32,2}

    function MC(params::Dict)
        Lx = params["Lx"]
        Ly = get(params, "Ly", Lx)
        T = params["T"]
        return new(T, zeros(Lx, Ly))
    end
end

function init!(mc::MC, data::LoadLeveller.MCData)
    mc.spin = rand!(data.rng, Bool, (mc.Lx, mc.Ly)) * 2 - 1
    return nothing
end

function sweep!(mc::MC, data::LoadLeveller.MCData)
    Lx = size(mc.spins, 1)
    for n = 1:length(mc.spins)
        i = Int32(rand!(data.rng) * length(mc.spins))
        x = 1 + i % Lx
        y = 1 + i / Lx

        ratio = exp(
            -2.0 / mc.T *
            mc.spins[x, y] *
            (mc.spins[x+1, y] + mc.spins[x-1, y] + mc.spins[x, y+1] + mc.spins[x, y-1]),
        )

        if ratio >= 1 || ratio > rand!(data.rng)
            mc.spins[1] *= -1
        end
    end
    return nothing
end

function write_checkpoint(mc::MC, d::HDF5.Group)
    d["spins"] = mc.spins
    return nothing
end

function read_checkpoint!(mc::MC, d::HDF5.Group)
    mc.spin = d["spins"]
    return nothing
end



end # module Ising
