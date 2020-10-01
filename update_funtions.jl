
function time_update(group,rooms,P)
    # counters to calculate incidence
    lat=0; pre=0; asymp=0; mild=0; hosp=0; sev=0; rec=0; dead=0;
    for x in group 
        x.tis += P.step 
        x.doi += P.step # increase day of infection. variable is garbage until person is latent

        if x.iso
            x.timeiso += P.step
            if x.timeiso >= max(P.iso_days,x.infp) ### a person must be recovered to get back 
                if x.health != DEAD 
                    if x.room_idx > 0
                        if rooms[x.room_idx].n_symp_res == 0 
                            x.iso = false
                        end
                    else
                        x.iso = false
                    end
                end
                #= if x.health == REC && x.room_idx > 0 ##if they are not dead and are an individual
                    rooms[x.room_idx].n_symp_res -= 1 ###This person is not consider symp anymore
                    #x.room_idx = x.reg_r_idx ##get back to their room
                end =#
            end
        end

        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
                :PRE  => begin move_to_pre(x); pre += 1; end
                :ASYMP => begin move_to_asymp(x); asymp += 1; end
                :MILD => begin move_to_mild(x,rooms); mild += 1; end
                :SEV  => begin move_to_sev(x,rooms); sev +=1; end  
                :HOSP  => begin move_to_hosp(x,rooms); hosp +=1; end    
                :REC  => begin move_to_recovered(x); rec += 1; end
                :DEAD  => begin move_to_dead(x); dead += 1; end
                _    => error("swap expired, but no swap set. $x")
            end

           
                
            #################################################################
            ##########Implement here the ones who got back from severe#######
            #################################################################
            #= if x.health == REC && x.wentTo == SEV
                x.iso = false
            elseif x.health == DEAD && x.room_idx > 0
                rooms[x.room_idx].n_av_beds -= 1 ###The bed is available
            end =#
        end
    end
    return lat, pre, asymp, mild, hosp, sev, rec, dead
end
export time_update


function move_to_latent(x::Humans)
    ## transfers human h to the incubation period and samples the duration
    x.health = LAT
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
    # the swap to asymptomatic is based on age group.
    # ask seyed for the references
    #asymp_pcts = (0.25, 0.25, 0.14, 0.07, 0.07)
    #symp_pcts = map(y->1-y,asymp_pcts) 
    #symp_pcts = (0.75, 0.75, 0.86, 0.93, 0.93) 
    
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.69, 0.71, 0.82, 0.87]
    age_thres = [4, 19, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)
     
    x.swap = rand() < (symp_pcts[g]) ? PRE : ASYMP 
    #x.got_inf = true
    ## in calibration mode, latent people never become infectious.
    if P.calibration 
        x.swap = LAT 
        x.exp = 999
    end
end
export move_to_latent

function move_to_asymp(x::Humans)
    ## transfers human h to the asymptomatic stage 
    x.health = ASYMP     
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
    x.swap = REC 
    x.wentTo = ASYMP
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp


function move_to_pre(x::Humans)
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period

    if x.room_idx < 0
        if rand() < θ[x.ag]
            x.swap = MILD
        else 
            x.swap = SEV
        end
    else
        x.swap = SEV
    end
    # calculate whether person is isolated
    #rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Humans,rooms::Array{Rooms,1})

    ## transfers human h to the mild infection stage for γ days
    x.health = MILD     
    x.tis = 0 
    x.exp = x.dur[4]
    
    x.swap = REC
    x.wentTo = MILD

    if !x.iso
        x.infp = x.dur[4]
        iso_ind(x)
        x.iso_symp = true
    end

    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
  
end
export move_to_mild

function move_to_sev(x::Humans,rooms::Array{Rooms,1})
    ## transfers human h to the mild infection stage for γ days
    x.health = SEV    
    x.tis = 0 
    x.exp = rand() ###we consider that it can take one day to move this person
   

    x.wentTo = SEV        

    if x.room_idx > 0
        age_thres = [60;69;79;89;100]
        a = findfirst(y-> y >= x.age,age_thres)
        ha = [0.182;0.151;0.134;0.093;0.071] #hospital admission probability
        if rand() < ha[a]
            x.swap = HOSP
        else
            mh = [0.084;0.188;0.241;0.279;0.354]
            if rand() < mh[a]
                x.exp = x.dur[4]  
                x.swap = DEAD
            else
                x.exp = x.dur[4]  
                x.swap = REC
            end
        end

       #=  mh = [0.084;0.188;0.241;0.279;0.354]
        if rand() < mh[a]
            x.exp = x.dur[4]  
            x.swap = DEAD
        end =#
        #set_new_room(x,rooms)
    else
        mh = [0.01/5, 0.01/5, 0.0135/3, 0.01225/1.5, 0.04/2]
        h = x.comorbidity == 1 ? 0.4 : 0.09
        x.health = SEV
        x.swap = UNDEF
        x.tis = 0 
        if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
            #x.exp = time_to_hospital
            x.exp = 0.0#time_to_hospital    
            x.swap = HOSP#rand() < c ? ICU : HOS        
        else ## no hospital for this lucky (but severe) individual 
            if rand() < mh[x.ag]
                x.exp = x.dur[4]  
                x.swap = DEAD
            else 
                x.exp = x.dur[4]  
                x.swap = REC
            end
        end

    end 


    if !x.iso
        x.infp = x.dur[4]
        iso_ind(x)
        x.iso_symp = true
    end

    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
   # if x.iso || rand() < p.fmild
   #     x.swap = MISO  
   #     x.exp = p.τmild
   # end
end
export move_to_sev


function move_to_hosp(x::Humans,rooms::Array{Rooms,1})
    ## transfers human h to the mild infection stage for γ days
    x.health = HOSP    
    x.tis = 0 
    x.exp = x.dur[4]-x.exp
    
    x.swap = REC


    x.wentTo = HOSP   
    #h = x.comorbidity == 1 ? 0.4 : 0.09
    c = x.comorbidity == 1 ? 0.33 : 0.25

    if x.room_idx > 0 #for residents
        age_thres = [60;69;79;89;100]
        a = findfirst(y-> y >= x.age,age_thres)
        mh = [0.084;0.188;0.241;0.279;0.354]
        x.exp = rand(truncated(Gamma(4.5, 2.75), 8, 17))
        if rand() < mh[a]
            x.swap = DEAD
        else
            x.swap = REC
        end
        rooms[x.room_idx].n_symp_res -= 1
    else #for HCW
        age_thres = [24;34;44;54;64;74;84;999]
        g = findfirst(y-> y >= x.age,age_thres)
        mh = [0.0005, 0.0022, 0.0057, 0.0160, 0.0401, 0.0696, 0.0893, 0.11]#non-icu service
        mc = [0.0009,0.0045,0.0115,0.0319,0.0801,0.1392,0.1786,0.22]#icu service +- 2x non-icu
        if rand() < c#[x.ag] ##goes to ICU
            x.exp = rand(truncated(Gamma(4.5, 2.75), 8, 17))+ 2
            if rand() < mc[x.ag]
                x.tis = x.exp  
                x.swap = DEAD
            end
        else
            x.exp = rand(truncated(Gamma(4.5, 2.75), 8, 17))
            if rand() < mh[x.ag]
                x.tis = x.exp
                x.swap = DEAD
            end
        end
    end
    x.infp = x.exp
  #=
    if rand() < mh[x.ag]
        x.tis = x.exp  
        x.swap = DEAD
    end =#



end
export move_to_hosp


function move_to_dead(x::Humans)
    # no level of alchemy will bring someone back to life. 
    x.health = DEAD
    x.swap = UNDEF
    x.tis = 0 
    x.exp = 999 ## stay recovered indefinitely
    x.iso_symp = false

    if (x.wentTo == SEV || x.h_t_delay > 0)  && x.room_idx > 0 #&& rooms[x.room_idx].n_symp_res > 0
        rooms[x.room_idx].n_symp_res -= 1
        #x.iso = false
    end
   # h.iso = true # a dead person is isolated
    #_set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(x::Humans)
    x.health = REC
    x.swap = UNDEF
    x.tis = 0
    #x.infp = 0 
    x.exp = 999 ## stay recovered indefinitely
    x.iso_symp = false

    if (x.wentTo == SEV || x.h_t_delay > 0)   && x.room_idx > 0 #&& rooms[x.room_idx].n_symp_res > 0
        
        rooms[x.room_idx].n_symp_res -= 1
    end

    #h.iso = false ## a recovered person has ability to meet others
    #_set_isolation(h, false)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function iso_ind(x::Humans)
    
   # println("$(x.idx) $(x.staff_type) $(x.h_t_delay) $(x.h_test)  $(x.health)")
    if P.iso_strat == :only
        x.iso = true
        x.timeiso = 0
        x.tested = false
        if x.room_idx > 0
            rooms[x.room_idx].n_symp_res += 1
        end
        x.iso_when = x.health
    elseif P.iso_strat == :total
        x.iso = true
        x.timeiso = 0
        x.tested = false
        if x.room_idx > 0
            rooms[x.room_idx].n_symp_res += 1
            
            pos = findall(y->y.room_idx == x.room_idx,residents)
            for kk = pos
                residents[kk].iso = true
                residents[kk].timeiso = 0
            end
        end
        x.iso_when = x.health
    end
    
end
function cleaning_rooms(rooms::Array{Rooms,1},pos)
    
    for i in pos
        rooms[i].contam = false
    end

end