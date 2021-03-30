using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.7       
    prov::Symbol = :ontario 
    calibration::Bool = false 
    modeltime::Int64 = 200
    initialinf::Int64 = 1
    fsevere::Float64 = 0.0 # fixed at 0.80
    fasymp::Float64 = 0.50 ## NOT USED ANYMORE ## percent going to asymp (may 10, removed all tests and references)
    fpre::Float64 = 1.0 ## NOT USED ANYMORE (percent going to presymptomatic)
    #frelasymp::Float64 = 0.11 ## relative transmission of asymptomatic
    #cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset
    number_of_residents::Int64 = 120
    hcw_per_shift::Int64 = 16
    staff_per_shift::Int64 = 5
    n_rooms::Int64 = 127
    test::Symbol = :np
    time_to_result::Int64 = 1

    #current_prev::Float64 = 0.0015

    n_shifts_pd::Int64 = 3 #shift_per_day
    n_hours_ps::Int64 = 8 #hours per shift

    number_of_bedrooms::Int64 = 0

    mild_red_idx::Float64 = 0.44
    asymp_red_idx::Float64 = 0.26 ####Changed here
    sev_red_idx::Float64 = 0.89

    step = 1/(n_shifts_pd*n_hours_ps)
    iso_days::Float64 = 14.0

    type_h::Symbol = :new
    iso_strat::Symbol = :total

    visitor_mask::Float64 = 0.3
    staff_mask::Float64 = 0.85
    nurse_mask::Float64 = 0.85

    normal_mask::Float64 = 0.67
    n95::Float64 = 0.95

    test_sens_hcw::Float64 = 1.0
    test_sens_res::Float64 = 1.0
    start_test::Int64 = 7
    test_interval::Int64 = 14
    testing_hcw::Bool = false
    testing_res::Bool = false

    #sub_hcw::Bool = false
    fixed_res::Int64 = 0

    vaccinating::Bool = false

    vac_cov_res::Float64 = 0.9
    vac_cov_hcw::Float64 = 0.7

    efficacy_red_min::Float64 = 0.0
    efficacy_red_max::Float64 = 0.0

    day_second_dose::Int64 = 14
    delay_immu::Int64 = 7

    prop_sev_res::Float64 = 0.8

    coming_inf::Bool = true
    prev_min::Float64 = 0.0005
    prev_max::Float64 = 0.001


    vac_efficacy_inf::Array{Array{Float64,1},1} = [[0.61],[0.61;0.935]]  #### 50:5:80
    vac_efficacy_symp::Array{Array{Float64,1},1} = [[0.921],[0.921;0.941]]  #### 50:5:80
    vac_efficacy_sev::Array{Array{Float64,1},1} = [[0.921],[0.921;1.0]]  #### 50:5:80

end

Base.@kwdef mutable struct ct_data_collect
    total_symp_id::Int64 = 0  # total symptomatic identified
    totaltrace::Int64 = 0     # total contacts traced
    totalisolated::Int64 = 0  # total number of people isolated
    iso_sus::Int64 = 0        # total susceptible isolated 
    iso_lat::Int64 = 0        # total latent isolated
    iso_asymp::Int64 = 0      # total asymp isolated
    iso_symp::Int64 = 0       # total symp (mild, inf) isolated
end

include("population.jl")
include("functions.jl")
include("update_funtions.jl")


## constants 
const residents = Array{Humans}(undef, 0)
const hcw = Array{Humans}(undef, 0)
const hcw_sub = Array{Humans}(undef, 0)
const rooms = Array{Rooms}(undef, 0) 
const P = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const ct_data = ct_data_collect()

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(P)
        setfield!(P, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

end

function main(P::ModelParameters,sim_idx::Int64)
    
    println(sim_idx)
    Random.seed!(256*sim_idx)
    reset_params(ip)

    creating_hosp_structure()
    
    aux = [rooms[i].n_beds for i = 1:length(rooms)]
    aux = cumsum(aux)

    if P.number_of_residents > 0
        resize!(residents,P.number_of_residents)
        for i = 1:P.number_of_residents
            ri = findfirst(x-> x >= i, aux)
            residents[i] = Humans()
            initializing_resident(residents[i],i,ri)
        end
    end
    
    
    hcw_per_shift,dist_PSW,dist_nurse,dist_HK,dist_diet = shift_numbers(P.type_h)
    
    creating_hcw_pop(hcw_per_shift,dist_PSW,dist_nurse,dist_HK,dist_diet)
    resize!(hcw_sub,length(hcw))
    for i = 1:length(hcw_sub)
        hcw_sub[i] = Humans()
      #=   for x in propertynames(hcw[i])
            setfield!(hcw_sub[i], x, getfield(hcw[i], x))
        end =#
    end

    if P.vaccinating
        vaccination_dose_1() 
    end

    first_inf = insert_infected(LAT, 1,hcw,rooms)

    time_length = P.modeltime*P.n_shifts_pd*P.n_hours_ps

    time_v = 1:time_length
    time_v = time_v*P.step

    lat_res_ct = zeros(Float64,time_length)
    lat_hcw_ct = zeros(Float64,time_length)
  
    asymp_res_ct = zeros(Float64,time_length)
    asymp_hcw_ct = zeros(Float64,time_length)
    
    pre_res_ct = zeros(Float64,time_length)
    pre_hcw_ct = zeros(Float64,time_length)
    
    mild_res_ct = zeros(Float64,time_length)
    mild_hcw_ct = zeros(Float64,time_length)
    
    hosp_res_ct = zeros(Float64,time_length)
    hosp_hcw_ct = zeros(Float64,time_length)
  
    sev_res_ct = zeros(Float64,time_length)
    sev_hcw_ct = zeros(Float64,time_length)
    
    rec_res_ct = zeros(Float64,time_length)
    rec_hcw_ct = zeros(Float64,time_length)
    
    dead_res_ct = zeros(Float64,time_length)
    dead_hcw_ct = zeros(Float64,time_length)
    

    t::Int64 = 1
    

    iso_tested_lat::Int64 = 0
    iso_tested_pre::Int64 = 0
    iso_tested_asymp::Int64 = 0

    total_lat::Int64 = 0
    total_pre::Int64 = 0
    total_asymp::Int64 = 0

    iso_tested_lat_res::Int64 = 0
    iso_tested_pre_res::Int64 = 0
    iso_tested_asymp_res::Int64 = 0


    iso_tested_lat = 0
    iso_tested_pre = 0
    iso_tested_asymp = 0 

    t_testing::Int64 = -P.start_test+1
    #t_testing = 0
    t_in_day::Int64 = 1
    #initiates the dynamics
    rnd::Float64 = (P.prev_max-P.prev_min)*rand()+P.prev_min
    for t_d = 1:P.modeltime ##run days
        if P.coming_inf
            if t_d%7 == 1
                rnd = (P.prev_max-P.prev_min)*rand()+P.prev_min
            end
            inserting_infections(rnd)
        end

        if P.vaccinating && t_d == (P.day_second_dose+P.delay_immu) 
            vaccination_dose_2_2()
        elseif P.vaccinating && t_d == (P.day_second_dose) 
            vaccination_dose_2()
        end

        if P.testing_res
            if t_d >= P.start_test
                #pos = findall(y->y.shift == n_shift,hcw)
                if t_testing%P.test_interval == 0
                    testing_individuals(residents)
                end
                lat_iso,pre_iso,asymp_iso = update_tested(residents)
                iso_tested_lat_res+=lat_iso
                iso_tested_pre_res+=pre_iso
                iso_tested_asymp_res+=asymp_iso
                
            end
        end
        #println(t_d)
        t_in_day = 1
            ##3number of contacts per day residents
        daily_contacts_res(residents)
       
        for n_shift = 1:P.n_shifts_pd #run the 3 shifts
            println(t_d,n_shift)
            #println("$t_d $n_shift")
            if P.testing_hcw
                if t_d >= P.start_test
                    pos = findall(y->y.shift == n_shift,hcw)
                    if t_testing%P.test_interval == 0
                        testing_individuals(hcw[pos])
                    end
                    lat_iso,pre_iso,asymp_iso = update_tested(hcw)
                    iso_tested_lat+=lat_iso
                    iso_tested_pre+=pre_iso
                    iso_tested_asymp+=asymp_iso
                    #t_testing += 1
                end
            end
            daily_contacts_hcw(n_shift)
            
            for h = 1:(P.n_hours_ps-1)#run the (n-1)th hours in a shift
               # println(t)
                contact_dynamics(residents,hcw,P,n_shift,t_in_day)
                
                (lat_res_ct[t], pre_res_ct[t], asymp_res_ct[t], mild_res_ct[t], hosp_res_ct[t], sev_res_ct[t], rec_res_ct[t], dead_res_ct[t],aux1,aux2,aux3) = time_update(residents,rooms,P) #updating the residents
                (lat_hcw_ct[t], pre_hcw_ct[t], asymp_hcw_ct[t], mild_hcw_ct[t], hosp_hcw_ct[t], sev_hcw_ct[t], rec_hcw_ct[t], dead_hcw_ct[t],aux1,aux2,aux3) = time_update(hcw,rooms,P) ##updating the hcw
                total_lat += aux1
                total_pre += aux2
                total_asymp += aux3
                (lat_iso, pre_iso, asymp_iso, mild_iso, hosp_iso, sev_iso, rec_iso, dead_iso,aux1,aux2,aux3) = time_update(hcw_sub,rooms,P) ##updating the hcw
                total_lat += aux1
                total_pre += aux2
                total_asymp += aux3
                
                lat_hcw_ct[t] += lat_iso
                pre_hcw_ct[t] += pre_iso
                asymp_hcw_ct[t] += asymp_iso
                mild_hcw_ct[t] += mild_iso
                hosp_hcw_ct[t] += hosp_iso
                sev_hcw_ct[t] += sev_iso
                rec_hcw_ct[t] += rec_iso
                dead_hcw_ct[t] += dead_iso
                t_in_day+=1
                t += 1
            end #end h
            #the 8-th hour is run here. 
            #Must copy everything inside the above loop here
            contact_dynamics(residents,hcw,P,n_shift,t_in_day)
            ## each 8 hours, we force one contact of residents with their roommates
            forcing_contact_res(residents,P)

            (lat_res_ct[t], pre_res_ct[t], asymp_res_ct[t], mild_res_ct[t], hosp_res_ct[t], sev_res_ct[t], rec_res_ct[t], dead_res_ct[t],aux1,aux2,aux3) = time_update(residents,rooms,P) #updating the residents
            (lat_hcw_ct[t], pre_hcw_ct[t], asymp_hcw_ct[t], mild_hcw_ct[t], hosp_hcw_ct[t], sev_hcw_ct[t], rec_hcw_ct[t], dead_hcw_ct[t],aux1,aux2,aux3) = time_update(hcw,rooms,P) ##updating the hcw
            total_lat += aux1
            total_pre += aux2
            total_asymp += aux3
            (lat_iso, pre_iso, asymp_iso, mild_iso, hosp_iso, sev_iso, rec_iso, dead_iso,aux1,aux2,aux3) = time_update(hcw_sub,rooms,P) ##updating the hcw
            total_lat += aux1
            total_pre += aux2
            total_asymp += aux3
            
            lat_hcw_ct[t] += lat_iso
            pre_hcw_ct[t] += pre_iso
            asymp_hcw_ct[t] += asymp_iso
            mild_hcw_ct[t] += mild_iso
            hosp_hcw_ct[t] += hosp_iso
            sev_hcw_ct[t] += sev_iso
            rec_hcw_ct[t] += rec_iso
            dead_hcw_ct[t] += dead_iso
            t_in_day+=1
            t += 1
            
        end #end n_shift
        t_testing += 1
    end #end t_d

    R0=length(findall(x->x.infected_by_type == hcw[first_inf[1]].staff_type,hcw))
    R0+=length(findall(x->x.infected_by_type == hcw[first_inf[1]].staff_type,residents))

    #return (lat_res_ct,lat_hcw_ct,pre_res_ct,pre_hcw_ct,asymp_res_ct,asymp_hcw_ct,mild_res_ct,mild_hcw_ct,sev_res_ct,sev_hcw_ct,hosp_res_ct,hosp_hcw_ct,rec_res_ct,rec_hcw_ct,dead_res_ct,dead_hcw_ct,R0,iso_tested_lat,iso_tested_pre,iso_tested_asymp,iso_tested_lat_res,iso_tested_pre_res,iso_tested_asymp_res)
    return (lat_res_ct,lat_hcw_ct,sum(pre_res_ct),sum(pre_hcw_ct),sum(asymp_res_ct),sum(asymp_hcw_ct),sum(sev_res_ct),sum(sev_hcw_ct),sum(hosp_res_ct),sum(hosp_hcw_ct),sum(dead_res_ct),sum(dead_hcw_ct),R0,iso_tested_lat,iso_tested_pre,iso_tested_asymp,total_lat,total_pre,total_asymp,sum(lat_res_ct),sum(lat_hcw_ct))

end
