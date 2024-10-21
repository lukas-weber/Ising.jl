# example_job.jl

using Carlo
using Carlo.JobTools
using Ising

tm = TaskMaker()

tm.sweeps = 20000
tm.thermalization = 2000
tm.binsize = 1

Ts = range(1, 3, 15)
tm.parallel_tempering = (mc = Ising.MC, parameter = :T, values = Ts, interval = 1)

Ls = [8, 12, 16]
for L in Ls
    tm.Lx = L
    tm.Ly = L
    task(tm)
end

job = JobInfo(
    splitext(@__FILE__)[1],
    ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = length(Ts),
)

start(job, ARGS)
