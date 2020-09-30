
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
