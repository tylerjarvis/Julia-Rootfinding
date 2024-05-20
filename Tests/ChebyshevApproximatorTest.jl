include("../../Julia-Rootfinding/yroots/ChebyshevApproximator.jl")

function transformPointsTest()
        # TESTING TRANSFORM
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

        @assert isapprox(transformPoints(x_1,a_1,b_1), result_1) "transform 1 dim failed"

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

        @assert isapprox(transformPoints(x_1,a_1,b_1), result_1) "transform 1 dim failed without array boxes"

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
                
        @assert isapprox(transformPoints(x_2,a_2,b_2), result_2) "transform 2 dim failed"


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

        @assert isapprox(transformPoints(x_3,a_3,b_3), result_3) "transform 3 dim failed"

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

        @assert isapprox(transformPoints(x_4,a_4,b_4), result_4) "transform 4 dim failed"
end

# RUNNING ALL TESTS
function allTests()
        transformPointsTest()
end

allTests()