include("../../Julia-Rootfinding/yroots/ChebyshevApproximator.jl")

function transformpoints_test()
    """Tests if the transformpoints functions
    correctly transforms points to the new bounds
    """
    # 1 dim test
    x_1 = [1;
            -1;
            .234775;
            -1/3]
    a_1 = [-5.54]
    b_1 = [1.2]
    result_1 = [1.2;
    -5.54;
    -1.37880825;
    -3.293333333333333]

    @assert isapprox(transformpoints(x_1,a_1,b_1), result_1) "transform 1 dim failed"

    # 1 dim test without arrays for upper/lower bounds
    x_1 = [1;
            -1;
            .234775;
            -1/3]
        a_1 = -5.54
    b_1 = 1.2
    result_1 = [1.2;
                -5.54;
                -1.37880825;
                -3.293333333333333]

    @assert isapprox(transformpoints(x_1,a_1,b_1), result_1) "transform 1 dim failed without array boxes"

    # 2 dim test
    x_2 = [1 1;
            -1 -1;
            .234775 .3456;
            -1/3 .99999999]
    a_2 = [-5.54 0]
    b_2 = [1.2 7.5]
    result_2 = [1.2 7.5;
            -5.54 0.0;
            -1.37880825 5.046;
            -3.293333333333333 7.4999999625]
            
    @assert isapprox(transformpoints(x_2,a_2,b_2), result_2) "transform 2 dim failed"


    # 3 dim test
    x_3 = [1 1 1;
            -1 -1 -1;
            -.2345555 .5 -.23455;
            -1/3 .99999999 1]
    a_3 = [0 0 1]
    b_3 = [1.2 7.5 4]
    result_3 = [1.2 7.5 4.;
            0. 0. 1.;
            0.4592667 5.625 2.148175;
            0.4 7.4999999625 4.] 

    @assert isapprox(transformpoints(x_3,a_3,b_3), result_3) "transform 3 dim failed"

    # 4 dim test
    x_4 = [1 1 1 1;
            -1 -1 -1 -1;
            4/9 -.9 5/166 -5/8;
            1/3 -9/10 1 -1]
    a_4 = [-19 4 7 -8]
    b_4 = [4/9 21/5 8 0]
    result_4 = [0.44444444444444 4.2 8 0;
            -19 4 7 -8;
            -4.956790123456791263834020355716 4.01 7.515060240963855164864071412012 -6.5;
            -6.037037037037038089692941866815 4.01 8 -8]

    @assert isapprox(transformpoints(x_4,a_4,b_4), result_4) "transform 4 dim failed"
end

# function hasconverged_test()

# end

function getfinal_degree_test()

end

function startedconverging_test()
    # .1 tolerance
    # should return true
    coefflist_1 = [3 2 1 .09999999999 .000456 .023 .01 .005]
    tol_1 = .1
    @assert startedconverging(coefflist_1,tol_1) == true ".1 tolerance that should have returned true that coefficients have started converging"

    # .1 tolerance
    # should return false
    coefflist_2 = [11 5.23 1.34 .09999999999 .000456 .023 .1 .005]
    tol_2 = .1
    @assert startedconverging(coefflist_2,tol_2) == false ".1 tolerance that should have returned false that coefficients have started converging"

    # .0000000000001 tolerance
    # should return true
    coefflist_3 = [75456 .0000000000000999978866666 .00000000000008726 .000000000000017382 .00000000000000327 .00000000000000019]
    tol_3 = .0000000000001
    @assert startedconverging(coefflist_3,tol_3) == true ".0000000000001 tolerance that should have returned true that coefficients have started converging"

    # .0000000000001 tolerance
    # should return false
    coefflist_4 = [75456 .000000000000100000000000003829 .00000000000008726 .000000000000017382 .00000000000000327 .00000000000000019]
    tol_4 = .0000000000001
    @assert startedconverging(coefflist_4,tol_4) == false ".0000000000001 tolerance that should have returned false that coefficients have started converging"
end

# RUNNING ALL TESTS
function alltests()
        transformpoints_test()
        # hasconverged_test()
        getfinal_degree_test()
        startedconverging_test()
end

alltests()