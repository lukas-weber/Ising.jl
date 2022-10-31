using Revise
using LoadLeveller
using Ising

LoadLeveller.run(Ising.MC, "test/job.data/parameters.json", LoadLeveller.SingleRunner)
