include("../../Julia-Rootfinding/yroots/ChebyshevApproximator.jl")

function test_transformpoints()
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

function test_getfinal_degree()
    # Test 1 checks if getFinalDegree properly returns a 5th degree 
    coeff1 = [4 8 .019 .09 .01 .001 .0001 .00005 .00004]
    tol1 = .02
    # Expected Outputs
    degree1 = 5
    epsval1 = .0002
    rho1 = 8.325532074018731520936853485182
    expected_outputs1 = [degree1 epsval1 rho1]
    degree1o, epsval1o, rho1o = getfinal_degree(coeff1,tol1)
    @assert isapprox(expected_outputs1, [degree1o epsval1o rho1o]) "getFinalDegree doesn't find 5th degree polynomial approximation correctly"


    # Test 2 checks if getFinalDegree properly finds 1 degree polynomials
    coeff2 = [.5 .02 .00111 .00287 .0028 .0015 .0005 .0004 .0015]
    tol2 = .01
    # Expected Outputs
    degree2 = 1
    epsval2 = .003
    rho2 = 12.909944487358055553727353981230
    expected_outputs2 = [degree2 epsval2 rho2]
    degree2o, epsval2o, rho2o = getfinal_degree(coeff2,tol2)
    @assert isapprox(expected_outputs2, [degree2o epsval2o rho2o]) "getFinalDegree doesn't find 1st degree polynomial approximation correctly"


    # Test 3 checks if getFinalDegree properly finds constant functions
    coeff3 = [.5 .001 .00111 .00287 .0028 .0015 .0005 .0004 .0015]
    tol3 = .01
    # Expected Outputs
    degree3 = 0
    epsval3 = .003
    rho3 = 166.666666666666657192763523198664
    expected_outputs3 = [degree3 epsval3 rho3]
    degree3o, epsval3o, rho3o = getfinal_degree(coeff3,tol3)
    @assert isapprox(expected_outputs3, [degree3o epsval3o rho3o]) "getFinalDegree doesn't return constant function appropriately"

    # Test 4 checks if degree of 1 is found when there are no non-zero coeff, but all coeff are still greater than the given tolerance
    coeff4 = [.0002 .0002 .0002 .0002 .0002 .0002 .0002 .0002 .0002]
    tol4 = .0001
    # # Expected Outputs
    degree4 = 1
    epsval4 = 0.0004
    rho4 = 0.707106781186547572737310929369
    expected_outputs4 = [degree4 epsval4 rho4]
    degree4o, epsval4o, rho4o = getfinal_degree(coeff4,tol4)
    @assert isapprox(expected_outputs4, [degree4o epsval4o rho4o]) "getFinalDegree doesn't find 1st degree polynomial approximation correctly"

    # Test 5 checks if getFinalDegree correctly handles the case where epsval will be 0 and the case when there are repeate maximum values in the coeff array that could be used in rho calculation
    coeff5 = [4 3 4 1 .1 .1 0 0 0]
    tol5 = .02
    # Expected Outputs
    degree5 = 5
    epsval5 = 4 * 1e-24
    rho5 = 9999.999999999994543031789362430573
    degree5o, epsval5o, rho5o = getfinal_degree(coeff5,tol5)
    @assert isapprox(epsval5o,epsval5) "getFinalDegree doesn't deal with epsval = 0 correctly"
    @assert isapprox(rho5o,rho5) "getFinalDegree doesn't deal with repeat max values correctly"
    @assert degree5 == degree5o "getFinalDegree doesn't find 5th degree correctly"

    # Test 6 verifies that getFinalDegree correctly evaluates a 6th degree polynomial when an array of length 11 is passed in.
    coeff6 = [4 8 6 3 .019 .09 .01 .001 .0001 .00005 .00004]
    tol6 = .02
    # Expected Outputs
    degree6 = 6
    epsval6 = .002
    rho6 = 3.98422018965844726424
    expected_outputs6 = [degree6 epsval6 rho6]
    degree6o, epsval6o, rho6o = getfinal_degree(coeff6,tol6)
    @assert isapprox(expected_outputs6, [degree6o epsval6o rho6o]) "getFinalDegree doesn't find 6th degree polynomial approximation correctly"
end

function test_startedconverging()
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

function test_check_constant_in_dimension()

    # Functions with bounds for testing
    f = (x_1,x_2,x_3,x_4) -> x_2 + 0.00000001*x_3
    bf = [5 5 1.1 5]
    af = [-5 -5 1 -5]
    
    g = (x_1,x_2,x_3,x_4,x_5) -> cos(x_1) + cos(x_2/80) + x_4 + sin(x_5)
    bg = [.5 .5 5 .5 5]
    ag = [0 0 -5 0 -5]

    # Tests with low tolerance
    tol = 1e-15
    @assert check_constant_in_dimension(f,af,bf,0,tol) == true "should return true since this dimension does not affect the function"
    @assert check_constant_in_dimension(f,af,bf,1,tol) == false "should return false since this dimension affects the function"
    @assert check_constant_in_dimension(f,af,bf,2,tol) == false "should return false since this dimension affects the funcion"
    @assert check_constant_in_dimension(f,af,bf,3,tol) == true "should return true since this dimension does not affect the function"
                                
    @assert check_constant_in_dimension(g,ag,bg,0,tol) == false "should return false since this dimension affects the function"
    @assert check_constant_in_dimension(g,ag,bg,1,tol) == false "should return false since this dimension affects the function"
    @assert check_constant_in_dimension(g,ag,bg,2,tol) == true "should return true since this dimension does not affect the function"
    @assert check_constant_in_dimension(g,ag,bg,3,tol) == false "should return false since this dimension affects the function"
    @assert check_constant_in_dimension(g,ag,bg,4,tol) == false "should return false since this dimension affects the function"

    # Tests with different tolerances where function changes little when selected dimension changes
    tol = 1e-10
    @assert check_constant_in_dimension(f,af,bf,2,tol) == false "should return false since this dimension affects the funcion"
    tol = 1e-7
    @assert check_constant_in_dimension(f,af,bf,2,tol) == true "should return true since the given change in this dimension is within the tolerance"

    tol = 1e-7
    @assert check_constant_in_dimension(g,ag,bg,1,tol) == false "should return true since the given change in this dimension is within the tolerance"
    tol = 1e-4
    @assert check_constant_in_dimension(g,ag,bg,1,tol) == true "should return false since this dimension affects the function"


    # Similar tolerance tests with the same functions, but less room for error
    tol = 1e-9
    @assert check_constant_in_dimension(f,af,bf,2,tol) == false "should return false since this dimension affects the funcion"
    tol = 1e-8
    @assert check_constant_in_dimension(f,af,bf,2,tol) == true "should return true since the given change in this dimension is within the tolerance"

    tol = 1e-6
    @assert check_constant_in_dimension(g,ag,bg,1,tol) == false "should return true since the given change in this dimension is within the tolerance"
    tol = 1e-5
    @assert check_constant_in_dimension(g,ag,bg,1,tol) == true "should return false since this dimension affects the function"

end

function test_has_converged()
    # Test case 1: Large tolerance, converged
    coeff = [0.1  0.0000001  -0.3]
    coeff2 = [0.11  -0.01  -0.31]
    tol = 0.01001
    @assert hasConverged(coeff, coeff2, tol) == true "has converged"

    # Test case 2: Small tolerance, not converged
    coeff = [0.1 0.2 0.3]
    coeff2 = [0.11 0.21 0.31]
    tol = 0.001
    @assert hasConverged(coeff, coeff2, tol) == false "has not converged"

    # Test case 3: Zero tolerance, not converged
    coeff = [0.1 0.2 0.3]
    coeff2 = [0.1000000001 0.2 0.3]
    tol = 0
    @assert hasConverged(coeff, coeff2, tol) == false "has not converged"

    # Test case 4: Large number of inputs, converged
    coeff = [i * 0.1 for i in 0:9]
    coeff2 = [(i * 0.1) + 0.0001 for i in 0:9]
    tol = 0.001
    @assert hasConverged(coeff, coeff2, tol) == true "has converged"
end

function test_create_meshgrid()
	#One dimensional array
	input = [5 1 2]
	expected = [5 1 2]
	got = create_meshgrid(input)
	@assert  got == expected "Test failed: One dimensional array \n\tExpected: $expected\n\tGot: $got"
	
	#One dimensional array with unnnecesarily nested values 
	input = [[[5]] [[1]] [[2]]]
	expected = [5 1 2]
	got = create_meshgrid(input)
	@assert  got == expected "Test failed: One dimensional array (nested) \n\tExpected: $expected\n\tGot: $got"

	#One dimensional array with unnnecesarily nested values 
	input = [[1 2] [3 4]]
	expected = [[1 1] [2 2] [3 4] [3 4]]
	got = create_meshgrid(input)
	@assert  got == expected "Test failed: 2x2 array \n\tExpected: $expected\n\tGot: $got"

	#Rectangular shape array 2 of length 3 
	input = [[1 2 3] [4 5 6]]
	expected = [[[1 1 1] [2 2 2] [3 3 3]]
				[[4 5 6] [4 5 6] [4 5 6]]]
	got = create_meshgrid(input)
	@assert  got == expected "Test failed: Rectangular array 3x2 \n\tExpected: $expected\n\tGot: $got"

	#Four dimensional array 
	input = [[1 2] [3 4] [5 6] [7 8]]
	expected = [[[[[1 1][1 1]][[1 1][1 1]]]
				[[[2 2][2 2]][[2 2][2 2]]]] 
				[[[[3 3][3 3]][[4 4][4 4]]]
				[[[3 3][3 3]][[4 4][4 4]]]] 
				[[[[5 5][6 6]][[5 5][6 6]]]
				[[[5 5][6 6]][[5 5][6 6]]]] 
				[[[[7 8][7 8]][[7 8][7 8]]]
				[[[7 8][7 8]][[7 8][7 8]]]]]
	got = create_meshgrid(input)
	@assert  got == expected "Test failed: four dimensional array \n\tExpected: $expected\n\tGot: $got"

end

function test_getApproxError()
	# Inputs
	degree1 = 5
	epsval1 = .0002
	rho1 = 8.325532074018731520936853485182

	degree2 = 1
	epsval2 = .003
	rho2 = 12.909944487358055553727353981230

	degree3 = 0
	epsval3 = .003
	rho3 = 166.666666666666657192763523198664

	degree4 = 1
	epsval4 = 0.0004
	rho4 = 0.707106781186547572737310929369

	degree5 = 5
	epsval5 = 4 * 1e-24
	rho5 = 9999.999999999994543031789362430573

	degree6 = 6
	epsval6 = .002
	rho6 = 3.98422018965844726424

	degs1 = [degree1]
	epsilons1 = [epsval1]
	rhos1 = [rho1]
	expected_approx_error1 = 0.00002730177111766866
	@assert isapprox(getApproxError(degs1,epsilons1,rhos1),expected_approx_error1) "wrong error value from 1 dim in getApproxError"

	degs2 = [degree1 degree2 degree3]
	epsilons2 = [epsval1 epsval2 epsval3]
	rhos2 = [rho1 rho2 rho3]
	expected_approx_error2 = 0.00183190901327267099
	@assert isapprox(getApproxError(degs2,epsilons2,rhos2),expected_approx_error2) "wrong error value getApproxError"

	degs3 = [degree1 degree2 degree3 degree4 degree5 degree6]
	epsilons3 = [epsval1 epsval2 epsval3 epsval4 epsval5 epsval6]
	rhos3 = [rho1 rho2 rho3 rho4 rho5 rho6]
	expected_approx_error3 = -0.87980728903851423972
	@assert isapprox(getApproxError(degs3,epsilons3,rhos3),expected_approx_error3) "wrong error value getApproxError"

	degs4 = [degree5]
	epsilons4 = [epsval5]
	rhos4 = [rho5]
	expected_approx_error4 = 4.0004000400040024e-28
	@assert isapprox(getApproxError(degs4,epsilons4,rhos4),expected_approx_error4) "wrong error value from getApproxError when error is very small"

end

# RUNNING ALL TESTS
function alltests()
    test_transformpoints()
    test_getfinal_degree()
    test_startedconverging()
    test_check_constant_in_dimension()
    test_has_converged()
    test_getApproxError()

    print("All tests finished")
end

