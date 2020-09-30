using Distributed
addprocs(4)

@everywhere using DelimitedFiles
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere using CSV
@everywhere include("parameters.jl")
@everywhere include("covid19_hosp_abm.jl")

function name_files(ip::ModelParameters)
    hospital_data_file = string("data/rooms_hosp.csv")
    return hospital_data_file
end


function create_folder()
    RF = string("results_$(ip.type_h)_$(ip.iso_strat)") ## 
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end

function dataprocess(results, ip::ModelParameters,NUMOFSIMS::Int64, fileappend="./")     
    ## takes the results of the pmap and stores it to file. 
    
    ## create empty vectors to store the results
    sim_time = ip.n_shifts_pd*P.n_hours_ps*P.modeltime
    resultsL_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsL_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsP_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsP_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsA_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsA_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsM_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsM_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsS_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsS_hcw = zeros(Int64, sim_time, NUMOFSIMS)


    resultsH_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsH_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsR_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsR_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    resultsD_res = zeros(Int64, sim_time, NUMOFSIMS)
    resultsD_hcw = zeros(Int64, sim_time, NUMOFSIMS)

    for i=1:NUMOFSIMS
        resultsL_res[:,i] = results[i][1]
        resultsL_hcw[:,i] = results[i][2]
        
        resultsP_res[:,i] = results[i][3]
        resultsP_hcw[:,i] = results[i][4]

        resultsA_res[:,i] = results[i][5]
        resultsA_hcw[:,i] = results[i][6]
        
        resultsM_res[:,i] = results[i][7]
        resultsM_hcw[:,i] = results[i][8]

        resultsS_res[:,i] = results[i][9]
        resultsS_hcw[:,i] = results[i][10]
        
        resultsH_res[:,i] = results[i][11]
        resultsH_hcw[:,i] = results[i][12]

        resultsR_res[:,i] = results[i][13]
        resultsR_hcw[:,i] = results[i][14]

        resultsD_res[:,i] = results[i][15]
        resultsD_hcw[:,i] = results[i][16]
        
    end
        
    writedlm(string("$fileappend", "_latent_res.dat"),  resultsL_res)
    writedlm(string("$fileappend", "_latent_hcw.dat"),  resultsL_hcw)
  
    writedlm(string("$fileappend", "_pre_res.dat"),  resultsP_res)
    writedlm(string("$fileappend", "_pre_hcw.dat"),  resultsP_hcw)
  
    writedlm(string("$fileappend", "_asymp_res.dat"),  resultsA_res)
    writedlm(string("$fileappend", "_asymp_hcw.dat"),  resultsA_hcw)
  
    writedlm(string("$fileappend", "_mild_res.dat"),  resultsM_res)
    writedlm(string("$fileappend", "_mild_hcw.dat"),  resultsM_hcw)
   
    writedlm(string("$fileappend", "_sev_res.dat"),  resultsS_res)
    writedlm(string("$fileappend", "_sev_hcw.dat"),  resultsS_hcw)
    

    writedlm(string("$fileappend", "_hosp_res.dat"),  resultsH_res)
    writedlm(string("$fileappend", "_hosp_hcw.dat"),  resultsH_hcw)

    writedlm(string("$fileappend", "_rec_res.dat"),  resultsR_res)
    writedlm(string("$fileappend", "_rec_hcw.dat"),  resultsR_hcw)
  
    writedlm(string("$fileappend", "_dead_res.dat"),  resultsD_res)
    writedlm(string("$fileappend", "_dead_hcw.dat"),  resultsD_hcw)
   
end


function runsim(simnum, ip::ModelParameters)
    
    folder = create_folder()
    dname = "$folder/beta_$(replace(string(ip.β), "." => "_"))"

    #rooms = creating_hosp_structure(hospital_data)
   # rooms1 = creating_hosp_structure_2(P) #creates the house structure that every simulation uses
    results = pmap(x-> main(ip,x),1:simnum)
    
    dataprocess(results,ip,simnum,dname)

end


#hospital_data = CSV.read("data/rooms_hosp.csv")

@everywhere ip = ModelParameters(β = 0.25, type_h = :new, iso_strat = :total) #new0.25 old 0.197
runsim(1000,ip)
