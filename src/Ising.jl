module Ising

using LoadLeveller
using HDF5
using Random

mutable struct MC <: LoadLeveller.AbstractMC
    T::Float64

    spins::Array{Int32,2}
end

function MC(params::Dict)
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)
    T = params[:T]
    return MC(T, zeros(Lx, Ly))
end

function LoadLeveller.init!(mc::MC, ctx::LoadLeveller.MCContext, params::AbstractDict)
    mc.spins = rand(ctx.rng, Bool, size(mc.spins)) * 2 .- 1
    return nothing
end

function periodic_elem(spins::AbstractArray, x::Integer, y::Integer)
    return spins[1+mod(x - 1, size(spins, 1)), 1+mod(y - 1, size(spins, 2))]
end

function LoadLeveller.sweep!(mc::MC, ctx::LoadLeveller.MCContext)
    Lx = size(mc.spins, 1)
    Ly = size(mc.spins, 2)


    for n = 1:length(mc.spins)
        i = floor(Int32, rand(ctx.rng) * length(mc.spins))
        x = 1 + i % Lx
        y = 1 + i ÷ Lx

        neighbor = (dx, dy) -> periodic_elem(mc.spins, x + dx, y + dy)
        ratio = exp(
            -2.0 / mc.T *
            mc.spins[x, y] *
            (neighbor(1, 0) + neighbor(-1, 0) + neighbor(0, 1) + neighbor(0, -1)),
        )

        if ratio >= 1 || ratio > rand(ctx.rng)
            mc.spins[x, y] *= -1
        end
    end
    return nothing
end

function LoadLeveller.measure!(mc::MC, ctx::LoadLeveller.MCContext)
    mag = sum(mc.spins) / length(mc.spins)

    energy = 0.0

    correlation = zeros(size(mc.spins, 1))
    for y = 1:size(mc.spins, 2)
        for x = 1:size(mc.spins, 1)
            neighbor = (dx, dy) -> periodic_elem(mc.spins, x + dx, y + dy)
            energy += -mc.spins[x, y] * (neighbor(1, 0) + neighbor(0, 1))
            # smart people use more lattice symmetries!
            correlation[x] += mc.spins[1, y] * mc.spins[x, y]
        end
    end

    LoadLeveller.measure!(ctx, :Magnetization, mag)
    LoadLeveller.measure!(ctx, :AbsMagnetization, abs(mag))
    LoadLeveller.measure!(ctx, :Magnetization2, mag^2)
    LoadLeveller.measure!(ctx, :Magnetization4, mag^4)

    LoadLeveller.measure!(ctx, :SpinCorrelation, correlation ./ size(mc.spins, 2))
    return nothing
end

function LoadLeveller.register_evaluables(
    ::Type{Ising.MC},
    eval::LoadLeveller.Evaluator,
    params::Dict,
)
    T = params[:T]
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)

    evaluate!(eval, :BinderRatio, [:Magnetization2, :Magnetization4]) do mag2, mag4
        return mag2 * mag2 / mag4
    end

    evaluate!(eval, :Susceptibility, [:Magnetization2]) do mag2
        return Lx * Ly * mag2 / T
    end

    evaluate!(eval, :SpinCorrelationK, [:SpinCorrelation]) do corr
        corrk = zero(corr)
        for i = 1:length(corr)
            for j = 1:length(corr)
                corrk[i] += corr[j] * cos(2 * pi / length(corr) * (i - 1) * (j - 1))
            end
        end
        return corrk
    end

    return nothing
end

function LoadLeveller.write_checkpoint(mc::MC, out::HDF5.Group)
    out["spins"] = mc.spins
    return nothing
end

function LoadLeveller.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.spins = read(in, "spins")
    return nothing
end



end # module Ising
