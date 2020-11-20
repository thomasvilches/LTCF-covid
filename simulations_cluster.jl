using Distributed
using ClusterManagers
#addprocs(4)
addprocs(SlurmManager(500), N=17, topology=:master_worker, exeflags="--project=.")

@everywhere using DelimitedFiles
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere using CSV
#@everywhere include("parameters.jl")
@everywhere include("covid19_hosp_abm.jl")

function name_files(ip::ModelParameters)
    hospital_data_file = string("data/rooms_hosp.csv")
    return hospital_data_file
end


function create_folder(ip::ModelParameters,end_file::String)
    #test = (ip.testing_res && ip.testing_hcw) ? "test_$(ip.test_interval)" : "not_testing"
    #RF = string("results_$(ip.type_h)_$(ip.iso_strat)_$(test)_$(ip.fixed_res)") ## 
    RF = string("/data/thomas-covid/LTCF/results_$(ip.type_h)_$(ip.iso_strat)_$end_file")
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

    resultsP_res = zeros(Int64, NUMOFSIMS)
    resultsP_hcw = zeros(Int64, NUMOFSIMS)

    resultsA_res = zeros(Int64, NUMOFSIMS)
    resultsA_hcw = zeros(Int64, NUMOFSIMS)

    resultsS_res = zeros(Int64, NUMOFSIMS)
    resultsS_hcw = zeros(Int64, NUMOFSIMS)

    resultsH_res = zeros(Int64, NUMOFSIMS)
    resultsH_hcw = zeros(Int64, NUMOFSIMS)

    resultsD_res = zeros(Int64, NUMOFSIMS)
    resultsD_hcw = zeros(Int64, NUMOFSIMS)

    results_R0 = zeros(Int64, NUMOFSIMS)
    resultsIso_lat = zeros(Int64, NUMOFSIMS)
    resultsIso_pre = zeros(Int64, NUMOFSIMS)
    resultsIso_asymp = zeros(Int64, NUMOFSIMS)
    resultsIso_lat_res = zeros(Int64, NUMOFSIMS)
    resultsIso_pre_res = zeros(Int64, NUMOFSIMS)
    resultsIso_asymp_res = zeros(Int64, NUMOFSIMS)

    resultsLT_res = zeros(Int64, NUMOFSIMS)
    resultsLT_hcw = zeros(Int64, NUMOFSIMS)
    for i=1:NUMOFSIMS
        resultsL_res[:,i] = results[i][1]
        resultsL_hcw[:,i] = results[i][2]
        
        resultsP_res[i] = results[i][3]
        resultsP_hcw[i] = results[i][4]

        resultsA_res[i] = results[i][5]
        resultsA_hcw[i] = results[i][6]
        
        resultsS_res[i] = results[i][7]
        resultsS_hcw[i] = results[i][8] 
        
        resultsH_res[i] = results[i][9]
        resultsH_hcw[i] = results[i][10]

        resultsD_res[i] = results[i][11]
        resultsD_hcw[i] = results[i][12]

        results_R0[i] = results[i][13]
        resultsIso_lat[i] = results[i][14]
        resultsIso_pre[i] = results[i][15]
        resultsIso_asymp[i] = results[i][16]

        resultsIso_lat_res[i] = results[i][17]
        resultsIso_pre_res[i] = results[i][18]
        resultsIso_asymp_res[i] = results[i][19]

        resultsLT_res[i] = results[i][20]
        resultsLT_hcw[i] = results[i][21]
        
    end
        
    writedlm(string("$fileappend", "_latent_res.dat"),  resultsL_res)
    writedlm(string("$fileappend", "_latent_hcw.dat"),  resultsL_hcw)
    
    writedlm(string("$fileappend", "_latent_t_res.dat"),  resultsLT_res)
    writedlm(string("$fileappend", "_latent_t_hcw.dat"),  resultsLT_hcw)
   
    writedlm(string("$fileappend", "_pre_res.dat"),  resultsP_res)
    writedlm(string("$fileappend", "_pre_hcw.dat"),  resultsP_hcw)
  
    writedlm(string("$fileappend", "_asymp_res.dat"),  resultsA_res)
    writedlm(string("$fileappend", "_asymp_hcw.dat"),  resultsA_hcw)
  
    #= writedlm(string("$fileappend", "_mild_res.dat"),  resultsM_res)
    writedlm(string("$fileappend", "_mild_hcw.dat"),  resultsM_hcw)
    =#
    writedlm(string("$fileappend", "_sev_res.dat"),  resultsS_res)
    writedlm(string("$fileappend", "_sev_hcw.dat"),  resultsS_hcw)
    

    writedlm(string("$fileappend", "_hosp_res.dat"),  resultsH_res)
    writedlm(string("$fileappend", "_hosp_hcw.dat"),  resultsH_hcw)

  #=    writedlm(string("$fileappend", "_rec_res.dat"),  resultsR_res)
    writedlm(string("$fileappend", "_rec_hcw.dat"),  resultsR_hcw) 
   =#
    writedlm(string("$fileappend", "_dead_res.dat"),  resultsD_res)
    writedlm(string("$fileappend", "_dead_hcw.dat"),  resultsD_hcw)

    writedlm(string("$fileappend", "_R0.dat"),  results_R0)
    writedlm(string("$fileappend", "_iso_asymp.dat"),  resultsIso_asymp)
    writedlm(string("$fileappend", "_iso_lat.dat"),  resultsIso_lat)
    writedlm(string("$fileappend", "_iso_pre.dat"),  resultsIso_pre)

    writedlm(string("$fileappend", "_iso_asymp_res.dat"),  resultsIso_asymp_res)
    writedlm(string("$fileappend", "_iso_lat_res.dat"),  resultsIso_lat_res)
    writedlm(string("$fileappend", "_iso_pre_res.dat"),  resultsIso_pre_res)
   
end


function runsim(simnum, ip::ModelParameters,end_file::String)
    
    println(end_file)
    folder = create_folder(ip,end_file)
    
    dname = "$folder/beta_$(replace(string(ip.β), "." => "_"))"

    #rooms = creating_hosp_structure(hospital_data)
   # rooms1 = creating_hosp_structure_2(P) #creates the house structure that every simulation uses
    results = pmap(x-> main(ip,x),1:simnum)
    
    dataprocess(results,ip,simnum,dname)

end


function run_all_scen(beta,tr1=[1;2;3],build=:new)

    @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = false,testing_res=false,fixed_res = 0,normal_mask = 0.0,n95 = 0.0,time_to_result = $tr) #new0.25 old 0.197
    runsim(8000,ip,"S0")
    
    @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = false,testing_res=false,fixed_res = 0,time_to_result = $tr) #new0.25 old 0.197
    runsim(8000,ip,"S1")
     
    @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = false,testing_res=false,fixed_res = 2,time_to_result = $tr) #new0.25 old 0.197
    runsim(8000,ip,"S2") 
    
    for tr in tr1
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:np,test_interval = 14,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3a-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:saliva,test_interval = 14,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3b-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:np,test_interval = 7,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3c-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:saliva,test_interval = 7,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3d-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:np,test_interval = 5, start_test = 5,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3e-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 0,test=:saliva,test_interval = 5, start_test = 5,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S3f-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:np,test_interval = 14,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4a-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:saliva,test_interval = 14,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4b-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:np,test_interval = 7,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4c-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:saliva,test_interval = 7,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4d-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:np,test_interval = 5, start_test = 5,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4e-$(ip.time_to_result)")
        
        @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = true,testing_res=false,fixed_res = 2,test=:saliva,test_interval = 5, start_test = 5,time_to_result = $tr) #new0.25 old 0.197
        runsim(8000,ip,"S4f-$(ip.time_to_result)")

    end
end 

function run_calibration(beta,tr=1,build=Symbol("new"))
    
    @everywhere ip = ModelParameters(β = $beta, type_h = :new, iso_strat = :total,testing_hcw = false,testing_res=false,fixed_res = 0,time_to_result = $tr) #new0.25 old 0.197
    runsim(8000,ip,"S1")
     
end
##################################################################33
############# old
###################################3

