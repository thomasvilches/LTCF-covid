


function matrix_p2p()
    #1 - normal residents
    #2- isolated residents
    #3 - hcw
    #4 - staff
    #=M = [
        8 0 10 3;
        0 0 8 4;
        10 8 8 5;
        3 4 5 10
    ]=#
    M = [
        8 0 10 0;
        0 0 8 0;
        10 8 8 0;
        0 0 0 0
    ]
    return M;
end


function matrix_p2r()
    #line
    #1 - normal residents
    #2 - isolated residents
    #3 - hcw
    #4 - staff

    #collumn
    #1 - normal rooms
    #2 - communal areas
    #3 - isolation rooms
    #4 - nursing station
    #=M = [
        2 3 0 1;
        0 3 0 1;
        3 2 5 2;
        5 5 3 4
    ]=#
    M = [
        2 3 0 1;
        0 3 0 1;
        3 2 5 2;
        0 0 0 0
    ]
    return M;
end



