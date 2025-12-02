#!/usr/bin/env julia

using Carlo
using Carlo.JobTools
using Ising

tm = TaskMaker()
tm.sweeps = 10000
tm.thermalization = 2000
tm.binsize = 5

tm.Lx = 20
tm.Ly = 20

tm.estimate_covariance = true # enable covariance!

tm.T = 2.27
task(tm)

job = JobInfo(splitext(@__FILE__)[1], Ising.MC;
    checkpoint_time="30:00",
    run_time="15:00",
    tasks=make_tasks(tm)
)

start(job, ARGS)

