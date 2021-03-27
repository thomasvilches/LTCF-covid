@enum HEALTH SUS LAT PRE ASYMP MILD SEV HOSP DEAD REC UNDEF
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]


Base.@kwdef mutable struct Humans
    idx::Int64 = 0 
    health::HEALTH = UNDEF
    swap::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    wentTo::HEALTH = UNDEF
    #nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Float64   = 0   # time in state 
    exp::Float64   = 999   # max statetime
    dur::NTuple{4, Float64} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Float64   = 999   # day of infection.
    #iso::Bool = false  ## isolated (limited contacts)
    room_idx::Int64 = -1
    reg_r_idx::Int64 = -1 ##this is an auxiliary int to store the room of the individual in case they go to isolation room, this was they can go back
    shift::Int64 = -1
    iso::Bool = false
    timeiso::Float64 = 0
    infp::Float64 = 0
    contacts_res::Array{Int64,1} = repeat([0],24)
    contacts_psw::Array{Int64,1} = repeat([0],24)
    contacts_nurse::Array{Int64,1} = repeat([0],24)
    contacts_hk::Array{Int64,1} = repeat([0],24)
    contacts_diet::Array{Int64,1} = repeat([0],24)
    infected_by::Int64 = -1
    infected_by_type::Symbol = :none
    comorbidity::Int16 = 0
    staff_type::Symbol = :none ##:psw :diet :hk :nurse
    mask_ef::Float64 = 0.0
    tested::Bool = false
    h_test::Int64 = -1
    h_t_delay::Int64 = -1
    iso_when::HEALTH = UNDEF
    tested_when::HEALTH = UNDEF
    iso_symp::Bool = false
    sub::Bool = false

    res_care::Array{Int64,1} = repeat([0],70)
    n_contacts::Int64 = 0
    contacts_done::Int64 = 0
    outside_inf::Bool = false

    vac_ef_symp::Float16 = 0.0 
    vac_ef_inf::Float16 = 0.0 
    vac_ef_sev::Float16 = 0.0
    vac_status::Int64 = 0
    red::Float64 = 0.0

end

Base.@kwdef mutable struct Rooms
    idx::Int64 = 0 
    typeofroom::String = "" 
    #typeofroom_backup::String = "" 
    n_beds::Int64 = -1
    n_av_beds::Int64 = -1
    contam::Bool = false
    n_symp_res::Int64 = 0
    ###Covid Hosp
end

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin        
        :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :canada => Distributions.Categorical(@SVector [0.0540, 0.1697, 0.3915, 0.2159, 0.1689])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        :ontario => Distributions.Categorical(@SVector [0.0519, 0.1727, 0.3930, 0.2150, 0.1674])
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])
        :newyork   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag