
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
    aux = cumsum(dist_HK)
    for k = 1:(sum(dist_HK))
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

        #= ###now for HCW we need to split the contacts
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
        end =#
    end
end

function daily_contacts_hcw(n_shift::Int64)
    
  
    aux_psw = findall(y-> y.staff_type == :psw && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_nurse = findall(y-> y.staff_type == :nurse && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_diet = findall(y-> y.staff_type == :diet && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_hk = findall(y-> y.staff_type == :hk && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    #aux_res = findall(y-> !(y.health in (DEAD,HOSP)),residents)

    for x in hcw
        x.contacts_done = 0
        x.n_contacts = 0
        for i = 1:length(x.contacts_res)
            x.contacts_res[i] = 0
            x.contacts_psw[i] = 0
            x.contacts_nurse[i] = 0
            x.contacts_hk[i] = 0
        end
        for i = 1:length(x.res_care)
            x.res_care[i] = 0
        end
    end

    psw_idx::Int64 = 0
    nurse_idx::Int64 = 0

    if P.fixed_res == 1
        aux = [i  for i = 1:length(residents)]
        for j in aux
            y = residents[j]
            psw_idx += 1
            psw_idx = psw_idx > length(aux_psw) ? 1 : psw_idx
            nurse_idx += 1
            nurse_idx = nurse_idx > length(aux_nurse) ? 1 : nurse_idx
            if !(y.health in (HOSP,DEAD))
                r = aux_psw[psw_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                
                r = aux_nurse[nurse_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                

            end
        end
    elseif P.fixed_res == 0
        aux = shuffle([i  for i = 1:length(residents)])

        for j in aux
            y = residents[j]
            psw_idx += 1
            psw_idx = psw_idx > length(aux_psw) ? 1 : psw_idx
            nurse_idx += 1
            nurse_idx = nurse_idx > length(aux_nurse) ? 1 : nurse_idx
            if !(y.health in (HOSP,DEAD))
                r = aux_psw[psw_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                
                r = aux_nurse[nurse_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                

            end
        end
    elseif P.fixed_res == 2

        aux = [i  for i = 1:length(residents)]

        aux_psw_1 = Int(floor(length(aux)/length(aux_psw)))
        aux_psw_2 = length(aux)%length(aux_psw)

        aux_nurse_1 = Int(floor(length(aux)/length(aux_nurse)))
        aux_nurse_2 = length(aux)%length(aux_nurse)

        n_c_psw = repeat([aux_psw_1],length(aux_psw))
        n_c_nurse = repeat([aux_nurse_1],length(aux_psw))
        for j = 1:aux_psw_2
            n_c_psw[j]+=1
        end
        for j = 1:aux_nurse_2
            n_c_nurse[j]+=1
        end
        

        psw_idx = 1
        nurse_idx = 1
        total = 0
        
        for j in aux
            y = residents[j]
            total += 1
            if total > cumsum(n_c_nurse)[nurse_idx]
               nurse_idx += 1
            end
            if total > cumsum(n_c_psw)[psw_idx]
                psw_idx += 1
            end
            
            if !(y.health in (HOSP,DEAD))
                r = aux_psw[psw_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                
                r = aux_nurse[nurse_idx]
                x = hcw[r]
                aux = x.shift
                hmin = (aux-1)*8+1
                hmax = aux*8
                
                r = rand(hmin:hmax)
                x.contacts_res[r] += 1
                x.n_contacts += 1
                x.res_care[x.n_contacts] = y.idx
                

            end
        end
    end

    
    
    for idx in aux_psw
        x = hcw[idx]
        aux = x.shift
        hmin = (aux-1)*8+1
        hmax = aux*8
       
        n = rand(2:4)
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
    for idx in aux_nurse
        x = hcw[idx]
        aux = x.shift
        hmin = (aux-1)*8+1
        hmax = aux*8
        
        n = rand(2:4)

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
    for idx in aux_diet
        x = hcw[idx]
        aux = x.shift
        hmin = (aux-1)*8+1
        hmax = aux*8
        
        n = 2
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
    for idx in aux_hk
        x = hcw[idx]
        aux = x.shift
        hmin = (aux-1)*8+1
        hmax = aux*8
        
        n = 2
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

function contact_dynamics(res::Array{Humans,1},hcw::Array{Humans,1},P::ModelParameters,c_shift::Int64,t::Int64)
    #contact matrices
    pos_res = findall(x -> !x.iso && x.health!=DEAD && x.health != HOSP,res)
    pos_iso = findall(x -> x.iso == true && x.health != DEAD && x.health != HOSP, res)

    pos_psw = findall(x -> x.shift == c_shift && x.staff_type == :psw && x.iso == false && x.health != DEAD && x.health != HOSP, hcw)
    pos_nurse = findall(x -> x.shift == c_shift && x.staff_type == :nurse && x.iso == false && x.health != DEAD && x.health != HOSP, hcw)
    pos_diet = findall(x -> x.shift == c_shift && x.staff_type == :diet && x.iso == false && x.health != DEAD && x.health != HOSP, hcw)
    pos_hk = findall(x -> x.shift == c_shift && x.staff_type == :hk && x.iso == false && x.health != DEAD && x.health != HOSP, hcw)
    
    
    totalinf_res::Int64 = 0
    totalinf_hcw::Int64 = 0
    totalinf_iso::Int64 = 0
    
    ih::HEALTH = UNDEF
    for i in pos_res
        x = residents[i]
        ih = x.health
        if x.health == ASYMP
           
            for j = 1:x.contacts_res[t]
                r = rand(pos_res)
                y = residents[r]
                if y.health == SUS
                    perform_contact(y,ih,P.asymp_red_idx,0.0,0.0)
                end
            end

        elseif x.health == PRE
            for j = 1:x.contacts_res[t]
                r = rand(pos_res)
                y = residents[r]
                if y.health == SUS
                    perform_contact(y,ih,0.0,0.0,0.0)
                end
            end
        end
    end
    red::Float64 = 0.0
    for i in [pos_nurse;pos_psw;pos_hk;pos_diet]
        x = hcw[i]
        ih = x.health
        if x.health == ASYMP
            red = P.asymp_red_idx
            for j = 1:x.contacts_nurse[t]
                r = rand(pos_nurse)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_psw[t]
                r = rand(pos_psw)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_diet[t]
                r = rand(pos_diet)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_hk[t]
                r = rand(pos_hk)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end
        elseif x.health == PRE
            red = 0.0
            for j = 1:x.contacts_nurse[t]
                r = rand(pos_nurse)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_psw[t]
                r = rand(pos_psw)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_diet[t]
                r = rand(pos_diet)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end

            for j = 1:x.contacts_hk[t]
                r = rand(pos_hk)
                y = hcw[r]
                if y.health == SUS
                    perform_contact(y,ih,red,P.normal_mask,P.normal_mask)
                end
            end
        end

    end #close for x

    for i in [pos_psw;pos_nurse]
        x = hcw[i]
        ih = x.health
       
        for kk = 1:x.contacts_res[t]

            x.contacts_done += 1
            #println("$(x.idx) $i $(x.contacts_done) $(x.res_care[x.contacts_done])")
            y = residents[x.res_care[x.contacts_done]]
            if x.health == SUS
                ih = y.health
                if y.health == ASYMP
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(x,ih,P.asymp_red_idx,red,0.0)
                elseif y.health == PRE
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(x,ih,0.0,red,0.0)

                elseif y.health == MILD
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(x,ih,P.mild_red_idx,red,0.0)
                elseif y.health == SEV
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(x,ih,P.sev_red_idx,red,0.0)
                end
            elseif x.health == ASYMP
                if y.health == SUS
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(y,ih,P.asymp_red_idx,red,0.0)
                end
            elseif x.health == PRE
                if y.health == SUS
                    red = y.iso ? P.n95 : P.normal_mask
                    perform_contact(y,ih,0.0,red,0.0)
                end
            end
        end
    end#close for x

end#close function

function perform_contact(sus_ind,ih,red_idx::Float64,mask_y::Float64,mask_x::Float64)

    r = rand()
    if r < P.β*(1-mask_y)*(1-mask_x)*(1-red_idx)
        sus_ind.swap = LAT
        sus_ind.exp = sus_ind.tis   ## force the move to latent in the next time step.
        sus_ind.sickfrom = ih ## stores the infector's status to the infectee's sickfrom
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

function testing_individuals(h::Array{Humans,1})

    for i = 1:length(h)
        x = h[i]
        
        ### I have the numbers for latent period t1, presymp, asymp, and symp
        if !x.iso_symp
            if x.health == LAT
                if x.swap == ASYMP
                    d = Int(floor(x.tis-x.dur[1]))
                    infec = infectivity(d)
                    if infec[1]<infec[2]
                        rnd = Distributions.Uniform(infec[1],infec[2])
                        r = rand(rnd)/2.0
                    else
                        r = 0
                    end
                    if rand() < r*P.test_sens
                        x.tested_when = x.health
                        x.tested = true
                        dist = [0.6;0.2;0.2]
                        a = findfirst(y-> rand() <= y,cumsum(dist))
                        x.h_t_delay =Int(a)
                    end
                else
                    d = Int(floor(x.tis-x.dur[1]-x.dur[3]))
                    infec = infectivity(d)
                    if infec[1]<infec[2]
                        rnd = Distributions.Uniform(infec[1],infec[2])
                        r = rand(rnd)
                    else
                        r = 0
                    end
                    
                    if rand() < r*P.test_sens
                        x.tested_when = x.health
                        x.tested = true
                        dist = [0.6;0.2;0.2]
                        a = findfirst(y-> rand() <= y,cumsum(dist))
                        x.h_t_delay =Int(a)
                        x.infp = x.dur[3]+x.dur[4]+x.dur[1]-x.tis
                    end

                end
            elseif  x.health == ASYMP
                d = Int(floor(x.tis))
                infec = infectivity(d)
                if infec[1]<infec[2]
                    rnd = Distributions.Uniform(infec[1],infec[2])
                    r = rand(rnd)/2.0
                else
                    r = 0
                end
                if rand() < r*P.test_sens
                    x.tested_when = x.health
                    x.tested = true
                    dist = [0.6;0.2;0.2]
                    a = findfirst(y-> rand() <= y,cumsum(dist))
                    x.h_t_delay = Int(a)
                end
            elseif x.health == PRE
                d = Int(floor(x.tis-x.dur[3]))
                infec = infectivity(d)
                if infec[1]<infec[2]
                    rnd = Distributions.Uniform(infec[1],infec[2])
                    r = rand(rnd)
                else
                    r = 0
                end
                
                if rand() < r*P.test_sens
                    x.tested_when = x.health
                    x.tested = true
                    dist = [0.6;0.2;0.2]
                    a = findfirst(y-> rand() <= y,cumsum(dist))
                    x.h_t_delay = Int(a)
                    x.infp = x.dur[3]+x.dur[4]-x.tis
                end
            end
        end
    end
end

function update_tested(h::Array{Humans,1})
    for x in h
        if x.tested && !x.iso_symp
            x.h_test += 1
            if x.h_test >= x.h_t_delay
                x.tested = false
                #x.h_test = -1
                #println("isolado por teste $(x.idx)")
                iso_ind(x)
                
            end
        end
    end
end


function infectivity(d::Int64)
    V = [-10;
    -9;
    -8;
    -7;
    -6;
    -5;
    -4;
    -3;
    -2;
    -1;
    0;
    1;
    2;
    3;
    4;
    5;
    6;
    7;
    8;
    9;
    10;
    11;
    12;
    13;
    14
    ]
    M =[0	0;
    0	0;
    0	0;
    0	0.02;
    0	0.03;
    0	0.09;
    0	0.22;
    0.05	0.43;
    0.13	0.67;
    0.58	1;
    0.92	1;
    0.69	0.95;
    0.46	0.8;
    0.29	0.62;
    0.18	0.49;
    0.11	0.37;
    0.06	0.24;
    0.03	0.12;
    0.02	0.09;
    0.01	0.05;
    0	0.03;
    0	0;
    0	0;
    0	0;
    0	0]

    if d >= minimum(V) && d <= maximum(V)
        a = findfirst(y->y==d,V)
    else
        a = 1
    end
 return M[a,:]
end


function create_subs(x::Humans)
   
    x.health = SUS
    x.swap = UNDEF
    x.sickfrom = UNDEF
    x.wentTo = UNDEF
    
    x.tis   = 0   # time in state 
    x.exp  = 999   # max statetime
    #x.doi::Float64   = 999   # day of infection.
    x.iso = false
    x.timeiso = 0
    x.iso_symp = false
    x.tested = false
    x.sub = true
    h_test = -1
    h_t_delay = -1
    iso_when = UNDEF
    tested_when = UNDEF
   
end