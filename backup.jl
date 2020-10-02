
function contact_dynamics(res::Array{Humans,1},rooms::Array{Rooms,1},hcw::Array{Humans,1},P::ModelParameters,c_shift::Int64)
    #contact matrices
    m_p2p = matrix_p2p()
    m_p2r = matrix_p2r()

    pos_res = findall(x->x.health in (SUS,LAT,PRE,ASYMP,REC),res)
    pos_hcw = findall(x -> x.shift == c_shift && x.health != MILD && x.health != SEV && x.health != DEAD,hcw)
    pos_iso = findall(x->x.health in (MILD,MISO),res)
    pos_staff = findall(x -> x.shift == c_shift && x.health != MILD && x.health != SEV && x.health != DEAD,staff)
    pos_rooms = findall(x->x.typeofroom=="bedroom",rooms) ##care here!!! parient/resident
    pos_com = findall(x->x.typeofroom=="communal",rooms) 
    pos_iso_r = findall(x->x.typeofroom=="isolation",rooms)
    pos_nst = findall(x->x.typeofroom=="nstation",rooms)



    totalinf_res::Int64 = 0
    totalinf_hcw::Int64 = 0
    totalinf_iso::Int64 = 0
    totalinf_staff::Int64 = 0
    #run all residents
    
    for i = 1:length(res)
        ### we need to make sure that every time a person leaves the hospital, the resident's health is set to UNDEF
        ##Running the three classes of individual contacts

        if res[i].health == SUS

            smc = sum(m_p2r[1,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2p[class,:]/smc)))

            #ROOMS
            totalinf_res += performing_contacts_r2p(res[i],rooms,pos_rooms,n_contacts_pg[1])
            #communal areas
            totalinf_res += performing_contacts_r2p(res[i],rooms,pos_com,n_contacts_pg[2])
            #Iso room (need it?)
            #totalinf_res += performing_contacts_r2p(ih,rooms,pos_iso_r,n_contacts_pg[3])

        elseif res[i].health in (ASYMP, PRE, MILD, MISO)
            class = res[i].iso ? 3 : 1
            smc = sum(m_p2p[class,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2p[class,:]/smc)))
            ih = res[i].health

            ##First group: normal Residents
            totalinf_res += performing_contacts_p2p(ih,res,pos_res,n_contacts_pg[1])

            ##second group: HCW
            totalinf_hcw += performing_contacts_p2p(ih,hcw,pos_hcw,n_contacts_pg[2])

            ##third group: isolated
            totalinf_iso += performing_contacts_p2p(ih,res,pos_iso,n_contacts_pg[3])

            ###Put here the environment contact

            smc = sum(m_p2r[class,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2r[class,:]/smc)))
            ih = res[i].health

            ##First group: Rooms
            performing_contacts_p2r(ih,rooms,pos_rooms,n_contacts_pg[1])

            ##second group: Communal areas
            performing_contacts_p2r(ih,rooms,pos_com,n_contacts_pg[2])
            ##third group: isolation (dont think we need this)
            #performing_contacts(ih,rooms,pos_iso_r,n_contacts_pg[3])
        end
        

    end
    
    ###find infectious for HCW
    
    for i = pos_hcw
        ### we need to make sure that every time a person leaves the hospital, the resident's health is set to UNDEF
        ##Running the three classes of individual contacts

        if hcw[i].health == SUS

            smc = sum(m_p2r[2,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2r[2,:]/smc)))

            #ROOMS
            totalinf_hcw += performing_contacts_r2p(hcw[i],rooms,pos_rooms,n_contacts_pg[1])
            #communal areas
            totalinf_hcw += performing_contacts_r2p(hcw[i],rooms,pos_com,n_contacts_pg[2])
            #Iso room (need it?)
            totalinf_hcw += performing_contacts_r2p(hcw[i],rooms,pos_iso_r,n_contacts_pg[3])

        elseif hcw[i].health in (ASYMP, PRE, MILD)
            smc = sum(m_p2p[2,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2p[2,:]/smc)))
            ih = hcw[i].health

            ##First group: normal Residents
            totalinf_res += performing_contacts_p2p(ih,res,pos_res,n_contacts_pg[1])

            ##second group: HCW
            totalinf_hcw += performing_contacts_p2p(ih,hcw,pos_hcw,n_contacts_pg[2])

            ##third group: isolated
            totalinf_iso += performing_contacts_p2p(ih,res,pos_iso,n_contacts_pg[3])

            ###Put here the environment contact

            smc = sum(m_p2r[2,:])
            n_contacts =  n_contacts_sample(smc,P)
            n_contacts_pg = Int.(round.(n_contacts*(m_p2r[2,:]/smc)))

            ##First group: Rooms
            performing_contacts_p2r(ih,rooms,pos_rooms,n_contacts_pg[1])
            ##second group: Communal areas
            performing_contacts_p2r(ih,rooms,pos_com,n_contacts_pg[2])
            ##third group: isolation (dont think we need this)
            performing_contacts_p2r(ih,rooms,pos_iso_r,n_contacts_pg[3])

        end
    end

end#close function









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

function daily_contacts_hcw(hcw::Array{Humans,1},n_shift::Int64)
    
    #dhr = [0.425000000, 0.241666667, 0.141666667, 0.050000000, 0.050000000, 0.050000000, 0.016666667, 0.016666667, 0.000000000, 0.008333333]
    #dhh = [0.50515464, 0.23711340, 0.00000000, 0.13402062, 0.00000000, 0.08247423, 0.00000000, 0.02061856, 0.00000000, 0.02061856]

    #dist_hr = Distributions.Categorical(dhr)
    #dist_hh = Distributions.Categorical(dhh)

    aux_psw = findall(y-> y.staff_type == :psw && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_nurse = findall(y-> y.staff_type == :nurse && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_diet = findall(y-> y.staff_type == :diet && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_hk = findall(y-> y.staff_type == :hk && y.shift == n_shift && (!y.iso || (y.iso && y.sub)),hcw)

    aux_res = findall(y-> !(y.health in (DEAD,HOSP)),residents)

    #= n_contacts_psw = [9;9;22]
    n_contacts_nurse = [32;32;64]
 =#
    N_psw = Int(ceil(length(aux_res)/length(aux_psw)))
    N_nurse = Int(ceil(length(aux_res)/length(aux_nurse)))
   
    for x in hcw
        
        for i = 1:length(x.contacts_res)
            x.contacts_res[i] = 0
            x.contacts_psw[i] = 0
            x.contacts_nurse[i] = 0
            x.contacts_hk[i] = 0
        end
    end
    
    
    for idx in aux_psw
        x = hcw[idx]
        aux = x.shift
        hmin = (aux-1)*8+1
        hmax = aux*8
        n_c::Int64 = N_psw

        for i = 1:n_c
            r = rand(hmin:hmax)
            x.contacts_res[r] += 1 
        end
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
        n_c::Int64 = N_psw

        for i = 1:n_c
            r = rand(hmin:hmax)
            x.contacts_res[r] += 1 
        end
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

function contact_dynamics(res::Array{Humans,1},rooms::Array{Rooms,1},hcw::Array{Humans,1},P::ModelParameters,c_shift::Int64,t_in_day::Int64)
    #contact matrices
    pos_res = findall(x -> !x.iso && x.health!=DEAD && x.health != HOSP,res)
    pos_iso = findall(x -> x.iso == true && x.health != DEAD && x.health != HOSP, res)

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


    for x in [pos_psw;pos_nurse]
        for i = 1:x.contacts_res[t]
            x.contacts_done += 1
            y = residents[x.res_care[x.contacts_done]]
            if x.health == SUS


            end

        end
    end

end#close function

function contact(x,res,hcw,pos_res::Array{Int64,1},pos_iso::Array{Int64,1},pos_psw::Array{Int64,1},pos_nurse::Array{Int64,1},pos_diet::Array{Int64,1},pos_hk::Array{Int64,1},line::Int64,t::Int64,P::ModelParameters)

    res_class = x.room_idx > 0 ? true : false
    
        
    if x.health == MILD || x.health == SEV
       
        reduction_idx = x.health == MILD ? P.mild_red_idx :  P.sev_red_idx
       
        performing_contacts_p2p(x,res,pos_res,x.contacts_res[t],P,reduction_idx,res_class,t) ##normal residents are in position 1
        #performing_contacts_p2p(x,res,pos_iso,x.contacts_res[t],P,reduction_idx,res_class,t) ##communal areas are in position 2
        performing_contacts_p2p(x,hcw,pos_psw,x.contacts_psw[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_nurse,x.contacts_nurse[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_diet,x.contacts_diet[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
        performing_contacts_p2p(x,hcw,pos_hk,x.contacts_hk[t],P,reduction_idx,!res_class,t) ##isolation rooms are in pos 3 of n contacts
       
    elseif x.health == ASYMP || x.health == PRE

        reduction_idx = x.health == ASYMP ? P.asymp_red_idx : 1.0
        
        if x.room_idx < 0 ####this contact part will be done in another function
            #performing_contacts_p2p(x,res,[pos_res;pos_iso],x.contacts_res[t],P,reduction_idx,res_class,t) ##normal residents are in position 1
        else
            performing_contacts_p2p(x,res,pos_res,x.contacts_res[t],P,reduction_idx,res_class,t) ##normal residents are in position 1
        end
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