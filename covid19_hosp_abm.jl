using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.7       
    prov::Symbol = :ontario 
    calibration::Bool = false 
    modeltime::Int64 = 300
    initialinf::Int64 = 1
    fsevere::Float64 = 0.0 # fixed at 0.80
    fasymp::Float64 = 0.50 ## NOT USED ANYMORE ## percent going to asymp (may 10, removed all tests and references)
    fpre::Float64 = 1.0 ## NOT USED ANYMORE (percent going to presymptomatic)
    frelasymp::Float64 = 0.11 ## relative transmission of asymptomatic
    #cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset
    number_of_residents::Int64 = 167
    hcw_per_shift::Int64 = 16
    staff_per_shift::Int64 = 5
    n_rooms::Int64 = 127

    n_shifts_pd::Int64 = 3 #shift_per_day
    n_hours_ps::Int64 = 8 #hours per shift

    number_of_bedrooms::Int64 = 0

    mild_red_idx::Float64 = 0.44
    asymp_red_idx::Float64 = 0.11
    sev_red_idx::Float64 = 0.89

    step = 1/(n_shifts_pd*n_hours_ps)
    iso_days::Float64 = 14.0

    type_h::Symbol = :new
    iso_strat::Symbol = :only

    visitor_mask::Float64 = 0.3
    staff_mask::Float64 = 0.85
    nurse_mask::Float64 = 0.85

    test_sens::Float64 = 1.0
    start_test::Int64 = 7
    test_interval::Int64 = 14
    testing::Bool = false

    sub_hcw::Bool = false

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

    insert_infected(LAT, 1,residents,rooms)

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
    

   # t::Int64 = 1
    t = 1
    
    #t_testing::Int64 = 0
    t_testing = 0
    #initiates the dynamics
    for t_d = 1:P.modeltime ##run days
        #t_in_day::Int64 = 1
        t_in_day = 1
            ##3number of contacts per day residents
        daily_contacts_res(residents)
       
        
        if P.testing
            if t_d >= P.start_test
                if t_testing%P.test_interval == 0
                    testing_individuals(hcw)
                end
                update_tested(hcw)
                t_testing += 1
            end
        end

        for n_shift = 1:P.n_shifts_pd #run the 3 shifts
            daily_contacts_hcw(hcw,n_shift)
            for h = 1:(P.n_hours_ps-1)#run the (n-1)th hours in a shift
               # println(t)
                contact_dynamics(residents,rooms,hcw,P,n_shift,t_in_day)
                (lat_res_ct[t], pre_res_ct[t], asymp_res_ct[t], mild_res_ct[t], hosp_res_ct[t], sev_res_ct[t], rec_res_ct[t], dead_res_ct[t]) = time_update(residents,rooms,P) #updating the residents
                (lat_hcw_ct[t], pre_hcw_ct[t], asymp_hcw_ct[t], mild_hcw_ct[t], hosp_hcw_ct[t], sev_hcw_ct[t], rec_hcw_ct[t], dead_hcw_ct[t]) = time_update(hcw,rooms,P) ##updating the hcw
                t_in_day+=1
                t += 1
            end #end h
            #the 8-th hour is run here. 
            #Must copy everything inside the above loop here
            contact_dynamics(residents,rooms,hcw,P,n_shift,t_in_day)
            ## each 8 hours, we force one contact of residents with their roommates
            forcing_contact_res(residents,P)

            (lat_res_ct[t], pre_res_ct[t], asymp_res_ct[t], mild_res_ct[t], hosp_res_ct[t], sev_res_ct[t], rec_res_ct[t], dead_res_ct[t]) = time_update(residents,rooms,P) #updating the residents
            (lat_hcw_ct[t], pre_hcw_ct[t], asymp_hcw_ct[t], mild_hcw_ct[t], hosp_hcw_ct[t], sev_hcw_ct[t], rec_hcw_ct[t], dead_hcw_ct[t]) = time_update(hcw,rooms,P) ##updating the hcw
            t_in_day+=1
            t += 1
        end #end n_shift

    end #end t_d

    return (lat_res_ct,lat_hcw_ct,pre_res_ct,pre_hcw_ct,asymp_res_ct,asymp_hcw_ct,mild_res_ct,mild_hcw_ct,sev_res_ct,sev_hcw_ct,hosp_res_ct,hosp_hcw_ct,rec_res_ct,rec_hcw_ct,dead_res_ct,dead_hcw_ct)
end