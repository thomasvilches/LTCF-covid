
function sample_epi_durations()
    # when a person is sick, samples the 
    
    pre_dist = truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = 
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods

    
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

function move_to_latent(x)
    ## transfers human h to the incubation period and samples the duration
    x.health = LAT
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 

    lat_dist = truncated(LogNormal(log(5.2), 0.1), 4, 7) # truncated between 4 and 7
    x.exp = rand(lat_dist) # get the latent period
    # the swap to asymptomatic is based on age group.
    # ask seyed for the references
    asymp_pcts = (0.25, 0.25, 0.14, 0.07, 0.07)    
    x.swap = rand() < asymp_pcts[x.ag] ? ASYMP : PRE 
    ## in calibration mode, latent people never become infectious.
    if P.calibration
        x.swap = LAT 
        x.exp = 999
    end
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = ASYMP     
    x.tis = 0 

    asy_dist = Gamma(5, 1)
    x.exp = Int.(ceil.(rand(asy_dist)))

    #x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period

    if rand() < θ[x.ag]
        x.swap = MILD
    else 
        x.swap = INF
    end
    # calculate whether person is isolated
    rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
    x.health = MILD     
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    if x.iso || rand() < p.fmild
        x.swap = MISO  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = MISO
    x.swap = REC
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    _set_isolation(x, true, :mi) 
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = INF
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    _set_isolation(x, false, :null) 
end
