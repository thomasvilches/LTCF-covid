using Plots
using Distributions
using FileIO
function hcwhcw()
    a = [2
    1
    1
    2
    1
    1
    1
    6
    3
    3
    5
    6
    4
    1
    1
    1
    1
    4
    3
    4
    4
    4
    1
    1
    1
    1
    2
    1
    1
    2
    1
    2
    2
    2
    2
    1
    1
    2
    2
    2
    3
    3
    1
    1
    1
    1
    2
    2
    2
    1
    4
    5
    2
    3
    3
    3
    4
    3
    1
    1
    3
    1
    3
    3
    2
    2
    2
    1
    2
    1
    2
    1
    2
    4
    2
    1
    1
    1
    1
    1
    1
    3
    1
    1
    1
    1
    1
    1
    1
    1
    2
    1
    1
    1
    1
    1
    1
    ]
    return a
end

a = hcwhcw()
maximum(a)
d = zeros(Float64,maximum(a))
map(x->d[x]+=1,a)

d = d/sum(d)
distr = Distributions.Categorical(d)
samp = rand(distr,10000)

b = Poisson(mean(a))
r = rand(b,10000)
for i in 1:length(r)
    if r[i] == 0
        r[i]+=1
    end
end


p=histogram(r,label = "poisson",norm = :probability,nbins = 0:1:15)
p=histogram!(samp,label = "data",norm = :probability,alpha=0.5,nbins = 0:1:15)
FileIO.save("hcwhc.png",p)


function hcwres()
    a= [
        4
        4
        2
        2
        2
        2
        1
        3
        2
        2
        7
        3
        4
        2
        3
        3
        2
        3
        2
        1
        1
        1
        1
        3
        2
        3
        2
        4
        3
        1
        2
        2
        1
        2
        2
        4
        1
        4
        1
        4
        2
        2
        4
        9
        6
        3
        2
        4
        7
        5
        9
        2
        1
        2
        7
        5
        6
        2
        4
        4
        11
        1
        7
        5
        3
        3
        3
        3
        5
        3
        2
        2
        2
        3
        2
        2
        1
        4
        2
        4
        2
        3
        5
        3
        3
        4
        3
        3
        7
        6
        8
        2
        3
        4
        1
        6
        2
        6
        3
        2
        4
        6
        1
        3
        4
        3
        3
        7
        5
        2
        3
        2
        3
        8
        2
        2
        3
        2
        2
        3        
    ]
    return a
end

a = hcwres()
d = zeros(Float64,maximum(a))
map(x->d[x]+=1,a)

d = d/sum(d)
distr = Distributions.Categorical(d)
samp = rand(distr,10000)

b = Poisson(mean(a))
r = rand(b,10000)
for i in 1:length(r)
    if r[i] == 0
        r[i]+=1
    end
end


p=histogram(r,label = "poisson",norm = :probability,nbins = 0:1:15)
p=histogram!(samp,label = "data",norm = :probability,alpha=0.5,nbins = 0:1:15)
FileIO.save("hcwres.png",p)

function resres()
    a = [1
    2
    2
    1
    1
    3
    2
    2
    1
    5
    8
    4
    5
    5
    5
    2
    4
    2
    1
    2
    3
    1
    1
    1
    1
    1
    2
    2
    5
    1
    3
    2
    2
    2
    3
    3
    2
    3
    3
    3
    3
    3
    3
    3
    3
    2
    3
    5
    5
    4
    3
    4
    5
    3
    4
    3
    5
    5
    3
    6
    9
    7
    8
    8
    6
    1
    4
    5
    5
    6
    8
    4
    2
    3
    1
    2
    4
    2
    2
    1
    2
    1
    5
    1
    1
    1
    1
    1
    1
    1
    4
    6
    6
    4
    4
    5
    2
    5
    3
    3
    5
    5]
    return a
end


a = resres()
d = zeros(Float64,maximum(a))
map(x->d[x]+=1,a)

d = d/sum(d)
distr = Distributions.Categorical(d)
samp = rand(distr,10000)

b = Poisson(mean(a))
r = rand(b,10000)
for i in 1:length(r)
    if r[i] == 0
        r[i]+=1
    end
end

p=histogram(r,label = "poisson",norm = :probability,nbins = 0:1:15)
p=histogram!(samp,label = "data",norm = :probability,alpha=0.5,nbins = 0:1:15)
FileIO.save("resres.png",p)