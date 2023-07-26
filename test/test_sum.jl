using Carlo.ResultTools
using Measurements
using StaticArrays
using DataFrames

function energy(spins)
    E = 0.0
    (Lx, Ly) = size(spins)
    for x = 1:Ly
        for y = 1:Lx
            E += -spins[x, y] * (spins[mod1(x + 1, Lx), y] + spins[x, mod1(y + 1, Ly)])
        end
    end
    return E
end

function solve(T::AbstractFloat, Lx::Integer, Ly::Integer; energy_cutoff = 10 * T)
    Z = 0.0
    Emean = 0.0
    M2 = 0.0
    M4 = 0.0

    E0 = -2 * Lx * Ly

    for spins in Iterators.product(((-1, 1) for _ = 1:Lx for _ = 1:Ly)...)
        E = energy(SMatrix{Lx,Ly}(spins))

        if E > energy_cutoff
            continue
        end

        w = exp(-(E - E0) / T)
        Z += w

        Emean += E * w
        m = float(sum(spins))
        M2 += m * m * w
        M4 += m * m * m * m * w
    end

    Emean /= Z * Lx * Ly
    M2 /= Z * (Lx * Ly)^2
    M4 /= Z * (Lx * Ly)^4

    return Dict(
        "Energy" => Emean,
        "Magnetization2" => M2,
        "Magnetization4" => M4,
        "BinderRatio" => M2^2 / M4,
    )

end


@testset "simple summation" begin
    tm = TaskMaker()
    tm.sweeps = 10000
    tm.thermalization = 2000
    tm.binsize = 500

    tm.Lx = 4
    tm.Ly = 4

    Ts = [0.1, 1.0, 3.0, 10.0]

    for T in Ts
        task(tm, T = T)
    end

    tasks = make_tasks(tm)

    mc_results = mktempdir() do dirname
        # dirname = "/tmp"
        job = JobInfo(
            dirname * "/test",
            Ising.MC;
            checkpoint_time = "15:00",
            run_time = "15:00",
            tasks = tasks,
        )

        start(Carlo.SingleScheduler, job)

        return DataFrame(ResultTools.dataframe(dirname * "/test.results.json"))
    end
    sum_results =
        DataFrame([solve(t.params[:T], t.params[:Lx], t.params[:Ly]) for t in tasks])

    obsnames = names(sum_results)
    for obsname in obsnames
        # hack for the case that the error is exactly zero
        mc_results_nudge =
            [m.val ± (m.err != 0 ? m.err : 1e-9) for m in mc_results[!, obsname]]
        z = stdscore.(mc_results_nudge, sum_results[!, obsname])
        chisqdof = sum(z .^ 2) / length(z)
        if chisqdof >= 5
            println(obsname)
            println(mc_results_nudge)
            println(sum_results[!, obsname])
        end
        @test chisqdof < 10
    end



end
