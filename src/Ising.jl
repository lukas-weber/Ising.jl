module Ising

using Carlo
using HDF5

mutable struct MC <: AbstractMC
    T::Float64

    spins::Array{Int32,2}
end

function MC(params::AbstractDict)
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)
    T = params[:T]
    return MC(T, zeros(Lx, Ly))
end

function Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
    mc.spins = rand(ctx.rng, Bool, size(mc.spins)) * 2 .- 1
    return nothing
end

function periodic_elem(spins::AbstractArray, x::Integer, y::Integer)
    return spins[mod1.((x, y), size(spins))...]
end

function Carlo.sweep!(mc::MC, ctx::MCContext)
    Lx = size(mc.spins, 1)

    for _ = 1:length(mc.spins)
        i = rand(ctx.rng, 0:length(mc.spins)-1)
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

function Carlo.measure!(mc::MC, ctx::MCContext)
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

    measure!(ctx, :Energy, energy / length(mc.spins))

    measure!(ctx, :Magnetization, mag)
    measure!(ctx, :AbsMagnetization, abs(mag))
    measure!(ctx, :Magnetization2, mag^2)
    measure!(ctx, :Magnetization4, mag^4)

    measure!(ctx, :SpinCorrelation, correlation ./ size(mc.spins, 2))
    return nothing
end

function Carlo.register_evaluables(::Type{Ising.MC}, eval::Evaluator, params::AbstractDict)
    T = params[:T]
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)

    evaluate!(eval, :BinderRatio, (:Magnetization2, :Magnetization4)) do mag2, mag4
        return mag2 * mag2 / mag4
    end

    evaluate!(eval, :Susceptibility, (:Magnetization2,)) do mag2
        return Lx * Ly * mag2 / T
    end

    evaluate!(eval, :SpinCorrelationK, (:SpinCorrelation,)) do corr
        corrk = zero(corr)
        for i = 1:length(corr), j = 1:length(corr)
            corrk[i] += corr[j] * cos(2π / length(corr) * (i - 1) * (j - 1))
        end
        return corrk
    end

    return nothing
end

function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["spins"] = mc.spins
    return nothing
end

function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.spins = read(in, "spins")
    return nothing
end


end # module Ising
