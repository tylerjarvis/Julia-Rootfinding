include("../../Julia-Rootfinding/src/ChebyshevApproximator.jl")
using Test

function test_all_ChebyshevApproximator()
    @testset "All tests in ChebyshevApproximatorTest.jl" begin
        test_transformPoints()
        test_getFinalDegree()
        test_startedConverging()
        test_checkConstantInDimension()
        test_hasConverged()
        test_createMeshgrid()
        test_getApproxError()
        test_intervalApproximateND()
        test_getChebyshevDegrees()
        test_chebApproximate()
    end
end

function test_transformPoints()
    @testset "transformPoints unit tests" begin
        # 1 dim test
        x_1 = [1;; -1;; 0.234775;; -1/3]
        a_1 = [-5.54]
        b_1 = [1.2]
        result_1 = [1.2;; -5.54;; -1.37880825;; -3.293333333333333]
        @test isapprox(transformPoints(x_1, a_1, b_1), result_1)
        # 1 dim test without arrays for upper/lower bounds
        x_1 = [1;; -1;; 0.234775;; -1/3]
        a_1 = -5.54
        b_1 = 1.2
        result_1 = [1.2;; -5.54;; -1.37880825;; -3.293333333333333]
        @test isapprox(transformPoints(x_1, a_1, b_1), result_1)
        # 2 dim test
        x_2 = [1;1;; -1;-1;; 0.234775;0.3456;; -1/3;0.99999999]
        a_2 = [-5.54;0]
        b_2 = [1.2;7.5]
        result_2 = [1.2;7.5;; -5.54;0.0;; -1.37880825;5.046;; -3.293333333333333;7.4999999625]
        @test isapprox(transformPoints(x_2, a_2, b_2), result_2)
        # 3 dim test
        x_3 = [1;1;1;; -1;-1;-1;; -0.2345555;0.5;-0.23455;; -1/3;0.99999999;1]
        a_3 = [0;0;1]
        b_3 = [1.2;7.5;4]
        result_3 = [1.2;7.5;4.;; 0.;0.;1.;; 0.4592667;5.625;2.148175;; 0.4;7.4999999625;4.]
        @test isapprox(transformPoints(x_3, a_3, b_3), result_3)
        # 4 dim test
        x_4 = [1;1;1;1;; -1;-1;-1;-1;; 4/9;-0.9;5/166;-5/8;; 1/3;-9/10;1;-1]
        a_4 = [-19;4;7;-8]
        b_4 = [4/9;21/5;8;0]
        result_4 = [0.44444444444444;4.2;8;0;; -19;4;7;-8;; -4.956790123456791263834020355716;4.01;7.515060240963855164864071412012;-6.5;; -6.037037037037038089692941866815;4.01;8;-8]
        @test isapprox(transformPoints(x_4, a_4, b_4), result_4)
    end
end

function test_getFinalDegree()
    @testset "getFinalDegree unit tests" begin
        # Test 1
        coeff1 = [4;8;0.019;0.09;0.01;0.001;0.0001;0.00005;0.00004;;]
        tol1 = 0.02
        degree1 = 5
        epsval1 = 0.0002
        rho1 = 8.325532074018731520936853485182
        expected_outputs1 = [degree1;epsval1;rho1]
        degree1o, epsval1o, rho1o = getFinalDegree(coeff1, tol1)
        @test isapprox(expected_outputs1, [degree1o;epsval1o;rho1o])
        # Test 2
        coeff2 = [0.5;0.02;0.00111;0.00287;0.0028;0.0015;0.0005;0.0004;0.0015]
        tol2 = 0.01
        degree2 = 1
        epsval2 = 0.003
        rho2 = 12.909944487358055553727353981230
        expected_outputs2 = [degree2;epsval2;rho2]
        degree2o, epsval2o, rho2o = getFinalDegree(coeff2, tol2)
        @test isapprox(expected_outputs2, [degree2o;epsval2o;rho2o])
        # Test 3
        coeff3 = [0.5;0.001;0.00111;0.00287;0.0028;0.0015;0.0005;0.0004;0.0015]
        tol3 = 0.01
        degree3 = 0
        epsval3 = 0.003
        rho3 = 166.666666666666657192763523198664
        expected_outputs3 = [degree3;epsval3;rho3]
        degree3o, epsval3o, rho3o = getFinalDegree(coeff3, tol3)
        @test isapprox(expected_outputs3, [degree3o;epsval3o;rho3o])
        # Test 4
        coeff4 = [0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002]
        tol4 = 0.0001
        degree4 = 1
        epsval4 = 0.0004
        rho4 = 0.707106781186547572737310929369
        expected_outputs4 = [degree4;epsval4;rho4]
        degree4o, epsval4o, rho4o = getFinalDegree(coeff4, tol4)
        @test isapprox(expected_outputs4, [degree4o;epsval4o;rho4o])
        # Test 5
        coeff5 = [4;3;4;1;0.1;0.1;0;0;0]
        tol5 = 0.02
        degree5 = 5
        epsval5 = 4.440892098500626e-16
        rho5 = 456.1401436878536
        degree5o, epsval5o, rho5o = getFinalDegree(coeff5, tol5)
        @test isapprox(epsval5o, epsval5)
        @test isapprox(rho5o, rho5)
        @test degree5 == degree5o
        # Test 6
        coeff6 = [4;8;6;3;0.019;0.09;0.01;0.001;0.0001;0.00005;0.00004]
        tol6 = 0.02
        degree6 = 6
        epsval6 = 0.002
        rho6 = 3.98422018965844726424
        expected_outputs6 = [degree6;epsval6;rho6]
        degree6o, epsval6o, rho6o = getFinalDegree(coeff6, tol6)
        @test isapprox(expected_outputs6, [degree6o;epsval6o;rho6o])
    end
end

function test_startedConverging()
    @testset "startedConverging unit tests" begin
        # Test 1
        coefflist_1 = [3;2;1;0.09999999999;0.000456;0.023;0.01;0.005]
        tol_1 = 0.1
        @test startedConverging(coefflist_1, tol_1) == true
        # Test 2
        coefflist_2 = [11;5.23;1.34;0.09999999999;0.000456;0.023;0.1;0.005]
        tol_2 = 0.1
        @test startedConverging(coefflist_2, tol_2) == false
        # Test 3
        coefflist_3 = [75456;0.0000000000000999978866666;0.00000000000008726;0.000000000000017382;0.00000000000000327;0.00000000000000019]
        tol_3 = 0.0000000000001
        @test startedConverging(coefflist_3, tol_3) == true 
        # Test 4
        coefflist_4 = [75456;0.000000000000100000000000003829;0.00000000000008726;0.000000000000017382;0.00000000000000327;0.00000000000000019]
        tol_4 = 0.0000000000001
        @test startedConverging(coefflist_4, tol_4) == false
    end
end

function test_checkConstantInDimension()
    @testset "checkConstantInDimension unit tests" begin
        # Functions with bounds for testing
        f = (x_1, x_2, x_3, x_4) -> x_2 + 0.00000001 * x_3
        bf = [5; 5; 1.1; 5]
        af = [-5; -5; 1; -5]
        g = (x_1, x_2, x_3, x_4, x_5) -> cos(x_1) + cos(x_2 / 80) + x_4 + sin(x_5)
        bg = [0.5; 0.5; 5; 0.5; 5]
        ag = [0; 0; -5; 0; -5]
        # Tests with low tolerance
        tol = 1e-15
        @test checkConstantInDimension(f, af, bf, 0, tol) == true
        @test checkConstantInDimension(f, af, bf, 1, tol) == false
        @test checkConstantInDimension(f, af, bf, 2, tol) == false
        @test checkConstantInDimension(f, af, bf, 3, tol) == true 
        @test checkConstantInDimension(g, ag, bg, 0, tol) == false 
        @test checkConstantInDimension(g, ag, bg, 1, tol) == false 
        @test checkConstantInDimension(g, ag, bg, 2, tol) == true
        @test checkConstantInDimension(g, ag, bg, 3, tol) == false 
        @test checkConstantInDimension(g, ag, bg, 4, tol) == false 
        # Tests with different tolerances where function changes little when selected dimension changes
        tol = 1e-10
        @test checkConstantInDimension(f, af, bf, 2, tol) == false
        tol = 1e-7
        @test checkConstantInDimension(f, af, bf, 2, tol) == true
        tol = 1e-7
        @test checkConstantInDimension(g, ag, bg, 1, tol) == false
        tol = 1e-4
        @test checkConstantInDimension(g, ag, bg, 1, tol) == true 
        # Similar tolerance tests with the same functions, but less room for error
        tol = 1e-9
        @test checkConstantInDimension(f, af, bf, 2, tol) == false
        tol = 1e-8
        @test checkConstantInDimension(f, af, bf, 2, tol) == true 
        tol = 1e-6
        @test checkConstantInDimension(g, ag, bg, 1, tol) == false 
        tol = 1e-5
        @test checkConstantInDimension(g, ag, bg, 1, tol) == true
    end
end

function test_hasConverged()
    @testset "hasConverged unit tests" begin
        # Test case 1
        coeff = [0.1; 0.0000001; -0.3]
        coeff2 = [0.11; -0.01; -0.31]
        tol = 0.01001
        @test hasConverged(coeff, coeff2, tol) == true
        # Test case 2
        coeff = [0.1; 0.2; 0.3]
        coeff2 = [0.11; 0.21; 0.31]
        tol = 0.001
        @test hasConverged(coeff, coeff2, tol) == false
        # Test case 3
        coeff = [0.1; 0.2; 0.3]
        coeff2 = [0.1000000001; 0.2; 0.3]
        tol = 0
        @test hasConverged(coeff, coeff2, tol) == false
        # Test case 4
        coeff = [i * 0.1 for i in 0:9]'
        coeff2 = [(i * 0.1) + 0.0001 for i in 0:9]'
        tol = 0.001
        @test hasConverged(coeff, coeff2, tol) == true
    end
end

function test_createMeshgrid()
    @testset "createMeshgrid unit tests" begin
        # One dimensional array
        input = [5; 1; 2]
        expected = [[5, 1, 2]]
        got = createMeshgrid(input)
        @test isapprox(got, expected)
        # 1x1x1 array 
        expected = [reshape([5], 1, 1, 1), reshape([1], (1, 1, 1)), reshape([2], (1, 1, 1))]
        got = createMeshgrid([5], [1], [2])
        @test got == expected 
        # One dimensional array with unnnecesarily nested values 
        expected = [[1; 1;; 2; 2], [3; 4;; 3; 4]]
        got = createMeshgrid([1; 2], [3; 4])
        @test got == expected 
        # Rectangular shape array 2 of length 3 
        expected = [[1; 1; 1;; 2; 2; 2;; 3; 3; 3], [4; 5; 6;; 4; 5; 6;; 4; 5; 6]]
        got = createMeshgrid([1; 2; 3], [4; 5; 6])
        @test got == expected 
        # Four dimensional array 
        expected = [[1; 1;; 1; 1;;; 1; 1;; 1; 1;;;; 2; 2;; 2; 2;;; 2; 2;; 2; 2], [3; 3;; 3; 3;;; 4; 4;; 4; 4;;;; 3; 3;; 3; 3;;; 4; 4;; 4; 4], [5; 5;; 6; 6;;; 5; 5;; 6; 6;;;; 5; 5;; 6; 6;;; 5; 5;; 6; 6], [7; 8;; 7; 8;;; 7; 8;; 7; 8;;;; 7; 8;; 7; 8;;; 7; 8;; 7; 8]]
        got = createMeshgrid([1; 2], [3; 4], [5; 6], [7; 8])
        @test got == expected 
    end
end

function test_getApproxError()
    @testset "getApproxError unit tests" begin
        # Inputs
        degree1 = 5
        epsval1 = 0.0002
        rho1 = 8.325532074018731520936853485182
        degree2 = 1
        epsval2 = 0.003
        rho2 = 12.909944487358055553727353981230
        degree3 = 0
        epsval3 = 0.003
        rho3 = 166.666666666666657192763523198664
        degree4 = 1
        epsval4 = 0.0004
        rho4 = 0.707106781186547572737310929369
        degree5 = 5
        epsval5 = 4 * 1e-24
        rho5 = 9999.999999999994543031789362430573
        degree6 = 6
        epsval6 = 0.002
        rho6 = 3.98422018965844726424
        degs1 = [degree1]
        epsilons1 = [epsval1]
        rhos1 = [rho1]
        expected_approx_error1 = 0.00002730177111766866
        @test isapprox(getApproxError(degs1, epsilons1, rhos1), expected_approx_error1)
        degs2 = [degree1; degree2; degree3]
        epsilons2 = [epsval1; epsval2; epsval3]
        rhos2 = [rho1; rho2; rho3]
        expected_approx_error2 = 0.00183190901327267099
        @test isapprox(getApproxError(degs2, epsilons2, rhos2), expected_approx_error2)
        degs3 = [degree1; degree2; degree3; degree4; degree5; degree6]
        epsilons3 = [epsval1; epsval2; epsval3; epsval4; epsval5; epsval6]
        rhos3 = [rho1; rho2; rho3; rho4; rho5; rho6]
        expected_approx_error3 = -0.87980728903851423972
        @test isapprox(getApproxError(degs3, epsilons3, rhos3), expected_approx_error3)
        degs4 = [degree5]
        epsilons4 = [epsval5]
        rhos4 = [rho5]
        expected_approx_error4 = 4.0004000400040024e-28
        @test isapprox(getApproxError(degs4, epsilons4, rhos4), expected_approx_error4)
    end
end

function test_intervalApproximateND()
    @testset "intervalApproximateND unit tests" begin
        f = (x, y, z) -> sin(5 * x + y + z)
        g = (x, y, z) -> cos(x * y * z)
        h = (x, y, z) -> x^2 + y^2 - z^2 - 1
        function1 = g
        degs1 =  [32; 5; 5]
        a1 =  [-4.2; 0.; 2.]
        b1 =  [3.3; 5.67; 3.3]
        expected_return_val1 = -0.215230482021798591452110827049
        expected_supnorm1 =  1.0
        return1, supnorm1 = intervalApproximateND(function1, degs1, a1, b1, true)
        @test isapprox(expected_return_val1, return1[2, 2, 3])
        @test isapprox(expected_supnorm1, supnorm1)
        function2 = f
        degs2 =  [32; 5; 5]
        a2 =  [-4.2; 0.; 2.]
        b2 =  [3.3; 5.67; 3.3]
        expected_return_val2 = 0.011960163308428820028161965183
        expected_supnorm2 = 0.999999291590031313958775172068
        return2, supnorm2 = intervalApproximateND(function2, degs2, a2, b2, true)
        @test isapprox(expected_return_val2, return2[2, 2, 3])
        @test isapprox(expected_supnorm2, supnorm2)
        function_3 = h
        degs_3 =  [2; 2; 17]
        a_3 =  [-4.2; 0.; 2.]
        b_3 =  [3.3; 5.67; 3.3]
        expected_return_val3 = 0
        expected_supnorm3 = 44.788899999999998158273228909820
        return3, supnorm3 = intervalApproximateND(function_3, degs_3, a_3, b_3, true)
        @test isapprox(expected_return_val3, return3[2, 2, 3]; atol = 2e-15)
        @test isapprox(expected_supnorm3, supnorm3)
        function4 = g
        degs4 =  [112; 75; 42]
        a4 =  [-4.2; 0.; 2.]
        b4 =  [3.3; 5.67; 3.3]
        expected_return_val4 = -0.012081908710314039068212110806
        return4 = intervalApproximateND(function4, degs4, a4, b4)
        @test isapprox(expected_return_val4, return4[2, 2, 3])
        # One dimensional tests
        oned_function = x -> x^3 + 3
        function5 = oned_function
        degs5 =  [8]
        a5 =  [-30.4]
        b5 =  [15.6]
        expected_return_val5 = 3041.75
        expected_supnorm5 = 28091.463999999996303813531994819641
        return5, supnorm5 = intervalApproximateND(function5, degs5, a5, b5, true)
        @test isapprox(expected_return_val5, return5[4])
        @test isapprox(expected_supnorm5, supnorm5)
        function6 = oned_function
        degs6 =  [3]
        a6 =  [-30.4]
        b6 =  [15.6]
        expected_return_val6 = 3041.75
        return6 = intervalApproximateND(function6, degs6, a6, b6, false)
        @test isapprox(expected_return_val6, return6[4])
        # Constant in a dimension test
        fx = (x, y) -> cos(x)
        function7 = fx
        degs7 =  [8; 0]
        a7 =  [-pi; -1e-5]
        b7 =  [pi; 1e-5]
        expected_return_val7 = -0.3042421775305829
        return5, supnorm5 = intervalApproximateND(function7, degs7, a7, b7, true)
        @test isapprox(expected_return_val7, return5[1])
    end
end

function test_getChebyshevDegrees()
    @testset "getChebyshevDegrees unit tests" begin
        relApproxTol = 1e-10
        f1 = x -> 4.2
        a1 = [-3.2]
        b1 = [4/3]
        expected_cheb_degs1 = [0]
        expected_epsilons1 = [0]
        expected_rhos1 = [Inf]
        cheb_degs1, epsilons1, rhos1 = getChebyshevDegrees(f1, a1, b1, relApproxTol)
        @test isapprox(expected_cheb_degs1, cheb_degs1)
        @test isapprox(expected_epsilons1, epsilons1)
        @test isapprox(expected_rhos1, rhos1)
        f2 = (x, y) -> 1e-5 * x + 1 + sin(y)
        a2 = [-1; -4.3]
        b2 = [1.2; 1]
        expected_cheb_degs2 = [1; 18]
        expected_epsilons2 = [0.00000000000000012262; 0.00000000000000017282]
        expected_rhos2 = [56243070.37339647859334945679; 6.55147385569233442482]
        cheb_degs2, epsilons2, rhos2 = getChebyshevDegrees(f2, a2, b2, relApproxTol)
        @test isapprox(expected_cheb_degs2, cheb_degs2)
        @test isapprox(expected_epsilons2, epsilons2)
        @test isapprox(expected_rhos2, rhos2)
        f3 = (x, y, z) -> 2 * cos(2 * x) + y * sin(y) + z
        a3 = [-10; -3; -4.3]
        b3 = [5; 2.3; 11/9]
        expected_cheb_degs3 = [41; 20; 1]
        expected_epsilons3 = [0.00000000000000007141; 0.00000000000000004932; 0.00000000000000006328]
        expected_rhos3 = [2.31102250026637268121; 5.49852805376172870666; 47683771.14014124125242233276]
        cheb_degs3, epsilons3, rhos3 = getChebyshevDegrees(f3, a3, b3, relApproxTol)
        @test isapprox(expected_cheb_degs3, cheb_degs3)
        @test isapprox(expected_epsilons3, epsilons3)
        @test isapprox(expected_rhos3, rhos3)
        f4 = (x1, x2, x3, x4) -> 1 + 7 * sin(1 / (13 * x2)) + 7 * x3 + x4
        a4 = [-10; 7e-5; -4.3; -2]
        b4 = [5; 2.3; 11/9; 1]
        expected_cheb_degs4 = [0; 74996; 1; 1]
        expected_epsilons4 = [0; 0.00000000000005330659; 0.00000000000000014466; 0.00000000000000028770]
        expected_rhos4 = [Inf; 1.00039596433691491129; 1855665305184964.25; 70038689.29537361860275268555]
        # cheb_degs4, epsilons4, rhos4 = getChebyshevDegrees(f4, a4, b4, relApproxTol)
        # @test isapprox(expected_cheb_degs4, cheb_degs4)
        # @test isapprox(expected_epsilons4, epsilons4)
        # @test isapprox(expected_rhos4, rhos4)
        # Same test as the previous, but this shouldn't throw a warning and the epsilon values should be larger since the relApproxTol passed in is larger 
        f5 = (x1, x2, x3, x4) -> 1 + 7 * sin(1 / (13 * x2)) + 7 * x3 + x4
        a5 = [-10; 7e-5; -4.3; -2]
        b5 = [5; 2.3; 11/9; 1]
        expected_cheb_degs5 = [0; 0; 1; 1]
        expected_epsilons5 = [0; 0.00080374052375353821; 0.00000000000000031672; 0.00000000000000057766]
        expected_rhos5 = [Inf; 519.95672738187363393081; 2542666157125790.00000000000000000000; 85612057.49236367642879486084]
        cheb_degs5, epsilons5, rhos5 = getChebyshevDegrees(f5, a5, b5, 1e-3)
        @test isapprox(expected_cheb_degs5, cheb_degs5)
        @test isapprox(expected_epsilons5, epsilons5)
        @test isapprox(expected_rhos5, rhos5)
    end
end

function test_chebApproximate()
    @testset "chebApproximate unit tests" begin
        f1 = (x,y,z) -> sin(5*x+y+z)
        a1 =  [-4.2; 0.; 2.]
        b1 =  [3.3; 5.67; 3.3]
        expected_error1 = 0.00000000000006049631
        expected_coeffs_shape1 = (23, 21, 48)
        expected_coeffs_val1a = 0.00000000000000000097
        expected_coeffs_val1b = -0.00000000000000000375

        coeffs1, error1 = chebApproximate(f1,a1,b1)
        @test isapprox(expected_error1,error1)
        @test (size(coeffs1) == expected_coeffs_shape1)
        @test isapprox(expected_coeffs_val1a,coeffs1[23,11,23])
        @test isapprox(expected_coeffs_val1b,coeffs1[23,21,48])

        f2 = x -> 3.1
        a2 = [-3.2]
        b2 = [4.2]
        expected_error2 = 0
        expected_coeffs_shape2 = (1,)
        expected_coeffs2 = [3.10000000000000008882]

        coeffs2, error2 = chebApproximate(f2,a2,b2)
        @test isapprox(expected_error2,error2)
        @test (size(coeffs2) == expected_coeffs_shape2)
        @test isapprox(expected_coeffs2,coeffs2)

        f3 = x -> 2*x^2-1
        a3 = [-4.3]
        b3 = [7/9]
        expected_error3 = 0
        expected_coeffs_shape3 = (3,)
        expected_coeffs3 = [11.64898148148147782877;-17.88506172839506191963;6.44595679012345623704]

        coeffs3,error3 = chebApproximate(f3,a3,b3)
        @test isapprox(expected_error3,error3)
        @test (size(coeffs3) == expected_coeffs_shape3)
        @test isapprox(expected_coeffs3,coeffs3)

        f4 = (x,y) -> 1e-5*x+1 + sin(y)
        a4 = [-1;-4.3]
        b4 = [1.2;1]
        expected_error4 = 0.00000000000000006226
        expected_coeffs_shape4 = (19, 2)
        expected_coeffs_val4 = -0.00000010641012795735

        coeffs4,error4 = chebApproximate(f4,a4,b4)
        @test isapprox(expected_error4,error4)
        @test (size(coeffs4) == expected_coeffs_shape4)
        @test isapprox(expected_coeffs_val4,coeffs4[13,1])

        f5 = (x1,x2,x3,x4) -> sin(7*x1/2) + x3^2 + x4
        a5 = [-1;3.1;-8/7;-13]
        b5 = [3.3;5.2;0;12]
        expected_error5 = 0.00000000000000073350
        expected_coeffs_shape5 = (2,3,1,29)
        expected_coeffs_val5a = 0
        expected_coeffs_val5b = -0.00000000000000013531

        coeffs5,error5 = chebApproximate(f5,a5,b5)
        @test isapprox(expected_error5,error5)
        @test (size(coeffs5) == expected_coeffs_shape5)
        @test isapprox(expected_coeffs_val5a,coeffs5[2,3,1,29])
        @test isapprox(expected_coeffs_val5b,coeffs5[1,3,1,21])

        f6 = (x1,x2,x3,x4) -> sin(7*x1/2) + x3^2 + x4
        a6 = [-1;3.1;-8/7;-13]
        b6 = [3.3;5.2;0;12]
        expected_error6 = 0.00000000001651317497
        expected_coeffs_shape6 = (2,3,1,24)
        expected_coeffs_val6a = -0.2133864467336392678
        expected_coeffs_val6b = -0.00000000000000077716

        coeffs6,error6 = chebApproximate(f6,a6,b6,1e-4)
        @test isapprox(expected_error6,error6)
        @test (size(coeffs6) == expected_coeffs_shape6)
        @test isapprox(expected_coeffs_val6a,coeffs6[1,1,1,1])
        @test isapprox(expected_coeffs_val6b,coeffs6[1,3,1,21])
    end
end
