In this version, each HCW takes care of a fixed number of individuals. The choosen individuals can be assigned acording to parameter fixed_res::Int64
fixed_res = 0 - randomly assignes the residents
fixed_res = 1 - the residents are fixed, but they are spread all over the house
fixed_res = 2 -  the residents are fixed and are taken from the same room, every time it is possible.

One can choose parameters for testing:
if staff will be tested; if residents will be tested; and sensitivity of test for both of them.
the interval of days between tests. And after how many days of introduction of infection the test starts.

File paper_plots.jl has some pre-programed analysis. one needs to choose beta and the type of building.

oct 27 - changes in probability of detection

Nov 28 - added a probability of HCW getting latent based on the current prevalence of cases in Canada (~0.001-0.005). Also, the results of the test is released after a fixed period that will be changed in order to check its influence in the system.

One can run "run_calibration(beta)" and use the file results_*_pre_res.dat to check the prevalence of symptomatic infections among residents and reach the reported prevalence. It is set for 120 days of epidemics.