WORK IN PROGRESS
This code is currently working.
In this version, each HCW takes care of a fixed number of individuals. The choosen individuals can be assigned acording to parameter fixed_res::Int64
fixed_res = 0 - randomly assignes the residents
fixed_res = 1 - the residents are fixed, but they are spread all over the house
fixed_res = 2 -  the residents are fixed and are taken from the same room, every time it is possible.
