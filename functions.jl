
function creating_hosp_structure()

     if P.type_h == :new
        dist = [57.1;42.8;0;0]
    elseif P.type_h == :old
        dist = [14;32;26;26]
    else
        println("error in building type")
        exit(0)
    end 

    dist = dist/sum(dist)
    aux = 1:4
    n_r = P.number_of_residents/sum(dist.*aux)
   
    dist = map(x-> Int(round(x)),dist*n_r)

    if sum(dist.*aux) < P.number_of_residents
        d = P.number_of_bedrooms - sum(dist.*aux)
        dist[1] += d
    end

    resize!(rooms, sum(dist))

    k::Int64 = 1
    for i = 1:length(dist)
        for j = 1:Int(dist[i])
            rooms[k] = Rooms()
            rooms[k].idx = k
            #rooms[i].typeofroom::String = "" 
            #rooms[i].#typeofroom_backup::String = "" 
            rooms[k].n_beds = i
            rooms[k].n_av_beds = 0
            rooms[k].contam = false
            rooms[k].n_symp_res = 0
            k+=1
        end
    end

    aux = 1:length(dist)
    n::Int64 = sum(dist.*aux)

    return n
    
end

function sample_epi_durations()
    # when a person is sick, samples the 
    lat_dist = Distributions.truncated(Gamma(3.122, 2.656),4,11.04) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

function initializing_resident(x,idx::Int64,ri::Int64)
    
    #agedist = get_province_ag(P.prov)
    x.health = SUS
    x.swap = UNDEF
    x.ag  = 5

    age_thres = [39;64;74;84;94;100]
    v = [0.0;6.6;11.4;27.3;43.9;10.8]

    v = v/sum(v)
    v = cumsum(v)
    v[1] = -0.01
    d = findfirst(y->rand()<= y,v)
    x.age = rand((age_thres[d-1]+1):age_thres[d])
    x.room_idx = ri
    #x.reg_r_idx = sr
    x.exp = 999
    #rooms[sr].n_av_beds -= 1
    x.idx = idx
    x.dur = sample_epi_durations()
    x.comorbidity = comorbidity(x.age)
end

function shift_numbers(type::Symbol)
    

    dist_PSW = [1/9.0;	1/9.0;	1/22.0]
    dist_nurse = [1/32.;1/32.0;	1/64.0]
    dist_diet = [1/32.;1/32.;1/32.]
    dist_HK = [1/32.;1/32.;1/32.]

    dist_PSW = map(x->Int(ceil(x)),P.number_of_residents*dist_PSW)
    dist_nurse = map(x->Int(ceil(x)),P.number_of_residents*dist_nurse)
    dist_diet = map(x->Int(ceil(x)),P.number_of_residents*dist_diet)
    dist_HK = map(x->Int(ceil(x)),P.number_of_residents*dist_HK)

    dist = dist_PSW+dist_nurse+dist_HK+dist_diet
    #= if type == :new
        #dist = [34;34;21]
       
    elseif type == :old
       
        dist = [72;61;21]
    else
        println("error in building type")
        exit(0)
    end =#

    return dist,dist_PSW,dist_nurse,dist_HK,dist_diet
end


function comorbidity(ag::Int16)
    agg = [4;17;39;59;79;999]
    g = findfirst(x->x>=ag,agg)
    prob = [0.05; 0.1;0.18; 0.38; 0.60; 0.72]
    com = rand() < prob[g] ? 1 : 0

    return com    
end

function initializing_hcw(x,P::ModelParameters,idx::Int64,si::Int64)
    
    agedist = get_province_ag(P.prov)
    hcw_adist = agedist.p[3:4]/sum(agedist.p[3:4])
    x.exp = 999
    x.health = SUS
    x.swap = UNDEF
    x.ag  = 2+findfirst(x->x>=rand(),cumsum(hcw_adist))
    x.age = rand(agebraks[x.ag]) 
    x.idx = idx
    x.shift = si#Int(ceil(idx/P.hcw_per_shift))
    x.dur = sample_epi_durations()
    x.comorbidity = comorbidity(x.age)
    
end

function creating_hcw_pop(hcw_per_shift::Array{Int64,1},dist_PSW::Array{Int64,1},dist_nurse::Array{Int64,1},dist_HK::Array{Int64,1},dist_diet::Array{Int64,1})
    
    nurse_ages = [20;24;34;54;64;75]
    nurse_dist_age = [6636;35566;67335;26203;5577]
    nurse_dist_age = Distributions.Categorical(nurse_dist_age/sum(nurse_dist_age))

    other_staff_age = [14;19;24;29;34;39;44;49;54;59;64;69;99]
    other_dist_age = [322.1;633.7;816.2;789.9;759.2;766.5;756.1;827.3;749.4;499.2;197.1;125.7]
    other_dist_age = Distributions.Categorical(other_dist_age/sum(other_dist_age))

    resize!(hcw,sum(hcw_per_shift))
    i::Int64 = 1
    ### Lets initialize PSW individuals
    aux = cumsum(dist_PSW)
    for k = 1:(sum(dist_PSW))
        hcw[i]  = Humans()
        ri = findfirst(x-> x >= k, aux)
        initializing_hcw(hcw[i],P,i,ri)
        hcw[i].staff_type = :psw
        r = rand(other_dist_age)
        hcw[i].age = rand(other_staff_age[r]+1:other_staff_age[r+1])
        hcw[i].ag = findfirst(y-> hcw[i].age in y,agebraks)
        hcw[i].mask_ef = P.staff_mask
        i += 1
    end
    #### Now, nurse
    aux = cumsum(dist_nurse)
    for k = 1:(sum(dist_nurse))
        hcw[i]  = Humans()
        ri = findfirst(x-> x >= k, aux)
        initializing_hcw(hcw[i],P,i,ri)
        hcw[i].staff_type = :nurse
        r = rand(nurse_dist_age)
        hcw[i].age = rand(nurse_ages[r]+1:nurse_ages[r+1])
        hcw[i].ag = findfirst(y-> hcw[i].age in y,agebraks)
        hcw[i].mask_ef = P.nurse_mask
        i += 1
    end
    ### Now dietary
    aux = cumsum(dist_diet)
    for k = 1:(sum(dist_diet))
        hcw[i]  = Humans()
        ri = findfirst(x-> x >= k, aux)
        initializing_hcw(hcw[i],P,i,ri)
        hcw[i].staff_type = :diet
        r = rand(other_dist_age)
        hcw[i].age = rand(other_staff_age[r]+1:other_staff_age[r+1])
        hcw[i].ag = findfirst(y-> hcw[i].age in y,agebraks)
        hcw[i].mask_ef = P.staff_mask
        i += 1
    end
    ### Now housekeeper
    aux = cumsum(dist_diet)
    for k = 1:(sum(dist_diet))
        hcw[i]  = Humans()
        ri = findfirst(x-> x >= k, aux)
        initializing_hcw(hcw[i],P,i,ri)
        hcw[i].staff_type = :hk
        r = rand(other_dist_age)
        hcw[i].age = rand(other_staff_age[r]+1:other_staff_age[r+1])
        hcw[i].ag = findfirst(y-> hcw[i].age in y,agebraks)
        hcw[i].mask_ef = P.staff_mask
        i += 1
    end


end

function insert_infected(health, num,res, rooms::Array{Rooms,1}) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS,res)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = res[i]
            if health == PRE 
                move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
            elseif health == LAT 
                move_to_latent(x)
            elseif health == MILD
                move_to_mild(x,rooms)
            elseif health == REC 
                move_to_recovered(x)
            else 
                error("can not insert human of health $(health)")
            end       
            x.sickfrom = UNDEF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
        end
    end    
    return h
end
export insert_infected

function daily_contacts_res(res::Array{Humans,1})
    
    drr = [0.411764706, 0.205882353, 0.107843137, 0.166666667, 0.049019608, 0.009803922, 0.039215686, 0.009803922]
    #drh = [0.43925234, 0.13084112, 0.06542056, 0.10280374, 0.07476636, 0.07476636, 0.05607477, 0.02803738, 0.02803738]

    dist_rr = Distributions.Categorical(drr)
    #dist_rh = Distributions.Categorical(drh)

    
    for x in res

        for i = 1:length(x.contacts_res)
            x.contacts_res[i] = 0
            x.contacts_psw[i] = 0
            x.contacts_nurse[i] = 0
            x.contacts_hk[i] = 0
        end
        n1 = Int(round(rand(dist_rr)/(0.36)))
        #n2 = Int(round(rand(dist_rh)/(0.63)))
        for i = 1:n1
            r = rand(1:24)
            x.contacts_res[r] += 1
        end

        ###now for HCW we need to split the contacts
        hmin = 1
        hmax = 8
        r = rand(hmin:hmax)
        x.contacts_psw[r] += 1
        r = rand(hmin:hmax)
        x.contacts_nurse[r] += 1
        hmin = 9
        hmax = 16
        r = rand(hmin:hmax)
        x.contacts_psw[r] += 1
        r = rand(hmin:hmax)
        x.contacts_nurse[r] += 1
        hmin = 17
        hmax = 24
        r = rand(hmin:hmax)
        x.contacts_psw[r] += 1
        r = rand(hmin:hmax)
        x.contacts_nurse[r] += 1

        if x.iso
            for i = 1:length(x.contacts_res)
                x.contacts_res[i] = 0#x.contacts_other_class[i]+x.contacts_same_class[i]
                #x.contacts_same_class[i] = 0
            end
        end
    end
end

function daily_contacts_hcw(hcw::Array{Humans,1})
    
    #dhr = [0.425000000, 0.241666667, 0.141666667, 0.050000000, 0.050000000, 0.050000000, 0.016666667, 0.016666667, 0.000000000, 0.008333333]
    #dhh = [0.50515464, 0.23711340, 0.00000000, 0.13402062, 0.00000000, 0.08247423, 0.00000000, 0.02061856, 0.00000000, 0.02061856]

    #dist_hr = Distributions.Categorical(dhr)
    #dist_hh = Distributions.Categorical(dhh)

    n_contacts_psw = [9;9;22]
    n_contacts_nurse = [32;32;64]
    for x in hcw
        
        for i = 1:length(x.contacts_res)
            x.contacts_res[i] = 0
            x.contacts_psw[i] = 0
            x.contacts_nurse[i] = 0
            x.contacts_hk[i] = 0
        end

        if !x.iso
            ###now for HCW we need to split the contacts
            aux = x.shift
            hmin = (aux-1)*8+1
            hmax = aux*8
            n_c::Int64 = 0

            if x.staff_type == :nurse
                n_c = n_contacts_nurse[x.shift]
            elseif x.staff_type == :psw
                n_c = n_contacts_psw[x.shift]
            end

            for i = 1:n_c
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1 
            end

            ### let's flip the copin and see the kind of hcw the hcw will meet

            if x.staff_type == :nurse || x.staff_type == :psw
                n = rand(2:4)
            else
                n = 2
            end

            aux = x.shift
            hmin = (aux-1)*8+1
            hmax = aux*8

            for j = 1:n
                r = rand(1:4)
                if r == 1
                    r = rand(hmin:hmax)
                    x.contacts_psw[r] += 1 
                    
                elseif r == 2
                    r = rand(hmin:hmax)
                    x.contacts_nurse[r] += 1 

                elseif r == 3
                    r = rand(hmin:hmax)
                    x.contacts_hk[r] += 1 
                else
                    r = rand(hmin:hmax)
                    x.contacts_diet[r] += 1 
                end
            end
        end
    end
end

function contact_dynamics(res::Array{Humans,1},rooms::Array{Rooms,1},hcw::Array{Humans,1},P::ModelParameters,c_shift::Int64,t_in_day::Int64)
    #contact matrices
    pos_res = findall(x -> !x.iso && x.health!=DEAD ,res)
    pos_iso = findall(x -> x.iso == true && (x.health != DEAD || x.health != HOSP), res)

    pos_psw = findall(x -> x.shift == c_shift && x.staff_type == :psw && x.iso == false && x.health != DEAD, hcw)
    pos_nurse = findall(x -> x.shift == c_shift && x.staff_type == :nurse && x.iso == false && x.health != DEAD, hcw)
    pos_diet = findall(x -> x.shift == c_shift && x.staff_type == :diet && x.iso == false && x.health != DEAD, hcw)
    pos_hk = findall(x -> x.shift == c_shift && x.staff_type == :hk && x.iso == false && x.health != DEAD, hcw)
    
    
    totalinf_res::Int64 = 0
    totalinf_hcw::Int64 = 0
    totalinf_iso::Int64 = 0
    
    pos_res_inf = findall(x -> !x.iso && x.health in (ASYMP,SEV,PRE,MILD) ,res)
    pos_iso_inf = findall(x -> x.iso == true && x.health in (ASYMP,SEV,PRE,MILD), res)
    
    pos_psw_inf = findall(x -> x.shift == c_shift && x.staff_type == :psw && x.iso == false && x.health in (ASYMP,SEV,PRE,MILD), hcw)
    pos_nurse_inf = findall(x -> x.shift == c_shift && x.staff_type == :nurse && x.iso == false && x.health in (ASYMP,SEV,PRE,MILD), hcw)
    pos_diet_inf = findall(x -> x.shift == c_shift && x.staff_type == :diet && x.iso == false && x.health in (ASYMP,SEV,PRE,MILD), hcw)
    pos_hk_inf = findall(x -> x.shift == c_shift && x.staff_type == :hk && x.iso == false && x.health in (ASYMP,SEV,PRE,MILD), hcw)

    #run all residents
    l_v = [length(pos_res_inf);length(pos_iso_inf);length(pos_psw_inf);length(pos_nurse_inf);length(pos_diet_inf);length(pos_hk_inf)] ##saving the length in a vector
    cl_v = cumsum(l_v) #saving the cummulative length
    n_total_ind = sum(l_v) ##number of updates that are performed

    for ww = 1:n_total_ind

        kk = rand(1:n_total_ind) #randomly takes one individual

        if kk <= cl_v[1] ###we can use the same test for non-iso and isolates
            i = pos_res_inf[kk]
            x = res[i]
            line = 1 ##contact matrix's line that we look at #not using anymore
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        elseif kk <= cl_v[2]
            kk -= cl_v[1]
            i = pos_iso_inf[kk]
            x = res[i]
            line = 2
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        elseif kk <= cl_v[3]
            kk -= (cl_v[2]) 
            i = pos_psw_inf[kk]
            x = hcw[i]
            line = 3
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        elseif kk <= cl_v[4]
            kk -= (cl_v[3]) 
            i = pos_nurse_inf[kk]
            x = hcw[i]
            line = 3
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        elseif kk <= cl_v[5]
            kk -= (cl_v[4]) 
            i = pos_diet_inf[kk]
            x = hcw[i]
            line = 3
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        elseif kk <= cl_v[6]
            kk -= (cl_v[5]) 
            i = pos_hk_inf[kk]
            x = hcw[i]
            line = 3
            contact(x,res,hcw,pos_res,pos_iso,pos_psw,pos_nurse,pos_diet,pos_hk,line,t_in_day,P)
        end
    end

end#close function

function contact(x,res,hcw,pos_res::Array{Int64,1},pos_iso::Array{Int64,1},pos_psw::Array{Int64,1},pos_nurse::Array{Int64,1},pos_diet::Array{Int64,1},pos_hk::Array{Int64,1},line::Int64,t::Int64,P::ModelParameters)

    res_class = x.room_idx > 0 ? true : false
    
    if x.health == SUS
        
    elseif x.health == MILD || x.health == SEV
       
        reduction_idx = x.health == MILD ? P.mild_red_idx :  P.sev_red_idx
        
       
        performing_contacts_p2p(x,res,pos_res,x.contacts_res[t],P,reduction_idx,res_class,t) ##normal residents are in position 1
        #performing_contacts_p2p(x,res,pos_iso,x.contacts_res[t],P,reduction_idx,res_class,t) ##communal areas are in position 2
        performing_contacts_p2p(x,hcw,pos_psw,x.contacts_psw[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_nurse,x.contacts_nurse[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_diet,x.contacts_diet[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_hk,x.contacts_hk[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
       
    elseif x.health == ASYMP || x.health == PRE

        reduction_idx = x.health == ASYMP ? P.asymp_red_idx : 1.0
        
        performing_contacts_p2p(x,res,pos_res,x.contacts_res[t],P,reduction_idx,res_class,t) ##normal residents are in position 1
        #performing_contacts_p2p(x,res,pos_iso,x.contacts_res[t],P,reduction_idx,res_class,t) ##communal areas are in position 2
        performing_contacts_p2p(x,hcw,pos_psw,x.contacts_psw[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_nurse,x.contacts_nurse[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_diet,x.contacts_diet[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_hk,x.contacts_hk[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
       
    end
end


function n_contacts_sample(average,P)
    dist = Poisson(average*P.step)
    nc = rand(dist)
    return nc
end


function performing_contacts_p2p(x::Humans,group,pos::Array{Int64,1},contacts::Int64,P::ModelParameters,reduction_idx::Float64,same_class::Bool,t::Int64)
    ih = x.health
    if length(pos) > 0
        for w = 1:contacts
            k = rand(pos)
            y = group[k]

            if x.room_idx > 0
                c = y.contacts_res[t]
                if c > 0 
                    y.contacts_res[t] -= 1
                end
            elseif x.staff_type == :nurse
                c = y.contacts_nurse[t]
                if c > 0 
                    y.contacts_nurse[t] -= 1
                end
            elseif x.staff_type == :psw
                c = y.contacts_psw[t]
                if c > 0 
                    y.contacts_psw[t] -= 1
                end
            elseif x.staff_type == :hk
                c = y.contacts_hk[t]
                if c > 0 
                    y.contacts_hk[t] -= 1
                end
            elseif x.staff_type == :diet
                c = y.contacts_diet[t]
                if c > 0 
                    y.contacts_diet[t] -= 1
                end
            end

            #c = same_class ? y.contacts_same_class[t] : y.contacts_other_class[t]
            if c > 0
               
                if y.health == SUS && y.swap == UNDEF
                    bf = P.β*reduction_idx ## baseline PRE
                    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
                    #=if ih == ASYMP
                        bf = bf * 0.11
                    elseif ih == MILD || ih == MISO 
                        bf = bf * 0.44
                    elseif ih == INF || ih == IISO 
                        bf = bf * 0.89
                    end=#
                    if rand() < bf*(1.0-x.mask_ef)*(1.0-y.mask_ef) ###implement  masking here
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = ih ## stores the infector's status to the infectee's sickfrom
                    end 
                end
            end
        end
    end
end

function forcing_contact_res(res::Array{Humans,1},P::ModelParameters)
    
    for i = 1:length(res)
        x = res[i]

        if x.health == MILD
            pos = findall(y->y.room_idx == x.room_idx, res)
            for j = pos
                y = res[j]
                if y.health == SUS && y.swap == UNDEF
                    bf = P.β*P.mild_red_idx
                    if rand() < bf
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = MILD ## stores the infector's status to the infectee's sickfrom
                    end 

                end
            end
        elseif x.health == ASYMP
            pos = findall(y->y.room_idx == x.room_idx, res)
            for j = pos
                y = res[j]
                if y.health == SUS && y.swap == UNDEF
                    bf = P.β*P.asymp_red_idx
                    if rand() < bf
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = ASYMP ## stores the infector's status to the infectee's sickfrom
                    end 

                end
            end
        elseif x.health == PRE
            pos = findall(y->y.room_idx == x.room_idx, res)
            for j = pos
                y = res[j]
                if y.health == SUS && y.swap == UNDEF
                    bf = P.β
                    if rand() < bf
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = PRE ## stores the infector's status to the infectee's sickfrom
                    end 

                end
            end
        elseif x.health == SEV
            pos = findall(y->y.room_idx == x.room_idx, res)
            for j = pos
                y = res[j]
                if y.health == SUS && y.swap == UNDEF
                    bf = P.β*P.sev_red_idx
                    if rand() < bf
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = SEV ## stores the infector's status to the infectee's sickfrom
                    end 

                end
            end

        end
    end

end