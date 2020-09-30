using DelimitedFiles
using Plots
using Parameters
using FileIO
include("parameters.jl") ###if you didnt run it before

ip = ModelParameters(β = 0.0058)

folder = "results_new_total/"

dname = "$folder/beta_$(replace(string(ip.β), "." => "_"))"
data_status_class = "_latent_res"

name_f = string(dname,data_status_class,".dat")

data = readdlm(name_f,header = false)

time_v = 1:size(data,1)
time_v = time_v*ip.step

pos = findall(x-> sum(data[:,x])>0,1:size(data,2))


m_v = map(x->mean(data[x,:]),1:size(data,1))

sum(m_v) ##area below the curve

#df = data[:,pos]
df = data
#m_v = map(x->mean(df[x,:]),1:size(df,1))
#Plots.plot(time_v,m_v,t = :bar,bins = 100,xlim = (1,10))

######################################################

df2 = zeros(Float64,ip.modeltime,size(df,2))

for i = 1:ip.modeltime
    for j = 1:size(df,2)
        base_l = ((i-1)*(ip.n_shifts_pd*ip.n_hours_ps)+1)
        base_h = ((i-1)*(ip.n_shifts_pd*ip.n_hours_ps)+ip.n_hours_ps*ip.n_shifts_pd)
        #for k = 1:(ip.n_shifts_pd*ip.n_hours_ps)
            df2[i,j] = sum(df[base_l:base_h,j])
        #end
    end
end


m_v = map(x->mean(df2[x,:]),1:size(df2,1))

sum(m_v) ##area below the curve

case = "total-new"

p1 = Plots.plot(m_v,t = :bar,xlab = data_status_class,label = "")
FileIO.save(string("barplot",data_status_class,"_",case,".png"),p1)
