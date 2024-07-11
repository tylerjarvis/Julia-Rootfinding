include("../../Julia-Rootfinding/src/ChebyshevSubdivisionSolver.jl")
using Test

function test_all_ChebyshevSubdivisionSolver()
    @testset "All tests in ChebyshevSubdivisionSolverTest.jl" begin
        test_getLinearTerms()
        test_linearCheck1()
        test_reduceSolveDim()
    end
end

function test_getLinearTerms()
    @testset "getLinearTerms unit tests" begin
        M_1 = reshape(collect(0:15),(2,2,2,2))
        A_1_expected = [8; 4; 2; 1]
        @test isapprox(A_1_expected,getLinearTerms(M_1))

        M_2 = reshape(collect(0:15),(2,2,4))
        A_2_expected = [4; 2; 1]
        @test isapprox(A_2_expected,getLinearTerms(M_2))

        M_3 = [43.2;12.2;-9.2]
        A_3_expected = [12.2]
        @test isapprox(A_3_expected,getLinearTerms(M_3))
        
        M_4 = reshape(collect(0:(2*15*4*6-1)),(6,4,15,2))
        A_4_expected = [360; 24; 6; 1]
        @test isapprox(A_4_expected,getLinearTerms(M_4))

        M_5 = reshape(collect(0:(3*4-1)),(4,1,3))
        A_5_expected = [4; 0; 1]
        @test isapprox(A_5_expected,getLinearTerms(M_5))


    end
end

function test_linearCheck1()
    @testset "linearCheck1 unit tests" begin
        A1 = [-1.12215095; -0.10896646; -0.20570858;;
                1.61496991;  0.36087537;  0.70263142;;
                1.11501566;  1.00547594; -1.96518304]
        totalErrs1 = [3.04770461; 7.78421744; 7.01361619]
        consts1 = [1.06289306; -0.36122354;  2.87185328]
        expecteda1 = [0.17843629682798045621; -5.97541921291522903914; 0.35379575125989282114]
        expectedb1 = [0.13891640768525181926; 0.26299355308293126399; 2.56893787868228251625]
        ao1, bo1 = linearCheck1(totalErrs1,A1,consts1)
        @test isapprox(expecteda1,ao1)
        @test isapprox(expectedb1,bo1)

        A2 = [1.80704361; 0.16184698; 0.12777497;;
                -1.50763366; -0.0496367; -0.03896966;;
                1.1471166; 1.60014553; -1.22405313]
        totalErrs2 = [3.12346641; 7.35040987; 6.71495483]
        consts2 = [0.43827528; 0.03844889; 2.66922695]
        expecteda2 = [-0.72849531284969915035; -3.19646507402361113037; -0.12454345016870327356]
        expectedb2 = [0.19996601042997674824; -0.13976516248493942030; 4.48583608458237481642]
        ao2, bo2 = linearCheck1(totalErrs2,A2,consts2)
        @test isapprox(expecteda2,ao2)
        @test isapprox(expectedb2,bo2)

        A3 = [1.80704361; 0.16184698; 0.12777497;;
                -1.50763366; -0.0496367;  -0.03896966;;
                -1.50763366; -0.0496367;  0.03896966]
        totalErrs3 = [3.12346641; 7.35040987; 6.71495483]
        consts3 = [0.43827528; 0.03844889; 2.66922695]
        expecteda3 = [0.08697917370722585417; -18.29888596005930878619; -23.44505688398909626358]
        expectedb3 = [0.24342093216001564615; 12.88296432840451899438; 16.58494523614445270709]
        ao3,bo3 = linearCheck1(totalErrs3,A3,consts3)
        @test isapprox(expecteda3,ao3)
        @test isapprox(expectedb3,bo3)

        A4 = reshape([1.50001956],(1,1)) #1 dimensional test
        totalErrs4 = [1.50001956]
        consts4 = [0.]
        expecteda4 = [0.]
        expectedb4 = [0.]
        ao4,bo4 = linearCheck1(totalErrs4,A4,consts4)
        @test isapprox(expecteda4,ao4)
        @test isapprox(expectedb4,bo4)

        A5 = [-.843750000; -1.87818750; -1.84275000; -.0107162044; -.00639690758; -.0137320288;;
                1.76080078; -.211296094; -.207309375; -.0159523129; .00889275657; -.0159523129;;
                .000176557928; -.00280220616; -.000701926758;  .0750781250; -1.28906250; -.0703125000;;
                1.36514364; -.163754473; -.160664766; -.163754473; .00504338969; -.00586371477;;
                .0289072432; -.240222558; -1.16359475; -1.12860484; -.190846439; -1.11307774;;
                .0903976590; -1.43472955; -.359386500;  0.0; 0.0;  0.0]
        totalErrs5 = [13.43636626; 31.17953809; 5.93851452; 33.20880558; 10.53793091; 1601.6394327]
        consts5 = [-0.07413729; -0.19634937;  0.90146117;-0.16101726; -0.72744493;  0.71811873]
        expecteda5 = [-14.92458223407407302830; -6.15390037469634965106; -6.29147538190204791420; -8.33713070909743869663; -2.20821696387878763090; -8.46738087673912076525]
        expectedb5 = [14.74884939851851761716; 6.07495480616285643549; 6.21101162935829531619; 7.04802595920109808958; 3.60684762763636346250; 7.16029350294975763802]
        ao5,bo5 = linearCheck1(totalErrs5,A5,consts5)
        @test isapprox(expecteda5,ao5)
        @test isapprox(expectedb5,bo5)

        A6 = [-1.09679818; -.939093750; -.921375000; -.0138146058; -.00106830799; -.0137307266;;
                .915137070; -1.05695980; -1.03701717; -.0149230031; .00237984138; -.0298460061;;
                .000607404435; -.00309112001; .000675877264;  .0750781250; -1.10873059; -.0379851973;;
                1.42079568; -1.63828770; -1.60737661; -1.63828770; -.00611697632; -.00852711940;;
                .901029633; -.240222558; -1.16359475; -1.12860484; -.0167700634; -1.01490072;;
                .310991071; -1.58265345;  .346049159;  0.0; 0.0; 0.0]
        totalErrs6 = [5.35450865; 14.31962505; 5.26765467; 30.10577724; 8.74980906; 1261.57287396]
        consts6 = [-0.06306681; -1.04453075; 1.44087973; -1.64179371; 0.23798175; 0.63996163]
        expecteda6 = [-3.88194523626944754113; -4.70178286246713916796; -4.81143253289919847759; -6.33103852363418972971; -1.15191610254029308535; -7.15236938643614461597]
        expectedb6 = [3.76694356841474720099; 4.56746866859671918348; 4.67453537376203964726; 6.75276584849662775412; 3.75106821937689982605; 7.62134481488987525211]
        ao6,bo6 = linearCheck1(totalErrs6,A6,consts6)
        @test isapprox(expecteda6,ao6)
        @test isapprox(expectedb6,bo6)

        A7 = [0; -.939093750; 0; -.0138146058; -.00106830799; -.0137307266;;
                0; -1.05695980; 0; -.0149230031; .00237984138; -.0298460061;;
                0; 0;  0;  0; 0; 0;;
                0; -1.63828770; 0; -1.63828770; -.00611697632; -.00852711940;;
                0; 0; 0; 0; 0; 0;;
                0; 0;  0;  0.0; 0.0; 0.0]
        totalErrs7 = [5.35450865; 14.31962505; 5.26765467; 30.10577724; 8.74980906; 1261.57287396]
        consts7 = [-0.06306681; -1.04453075; 1.44087973; -1.64179371; 0.23798175; 0.63996163]
        expecteda7 = [-Inf; -4.70178286246713916796; -Inf; -17.37636774053788002448; -4920.67627681775957171340; -388.96542615596172254300]
        expectedb7 = [Inf; 4.56746866859671918348; Inf; 15.37208764980656283683; 4383.87717081762366433395; 379.77919561809636661565]
        ao7,bo7 = linearCheck1(totalErrs7,A7,consts7)
        @test isapprox(expecteda7,ao7)
        @test isapprox(expectedb7,bo7)
    end
end

function test_reduceSolveDim()
    # 2d case
    Ms_1 = [[3.1;-2;34;.000001;.0002;;4;1;2;.0002;.0102], [9.836;3.1;-.031;.001;-.045;.00002;;15.2344;42.234;2;4;3;.0023]]
    errors_1 = [.00002;.0001]
    trackedInterval_1 = TrackedInterval([1;1;;-3;2])
    trackedInterval_1.reducedDims = [1]
    trackedInterval_1.solvedVals = [1.;2.;3.]
    dim_1 = 0
    expected_new_Ms_1 = [[3.10000000000000008882; -2.00000000000000000000; 34.00000000000000000000; 0.00000100000000000000; 0.00020000000000000001]]
    expected_new_errors_1 = [0.00002]
    expected_new_trackedInterval_1 = [-3;2]

    new_Ms_1, new_errors_1, new_trackedInterval_1 = reduceSolvedDim(Ms_1,errors_1,trackedInterval_1,dim_1)
    @test isapprox(expected_new_Ms_1[1],new_Ms_1[1])
    @test isapprox(expected_new_errors_1,new_errors_1)
    @test isapprox(expected_new_trackedInterval_1,new_trackedInterval_1.interval)
    @test (1 == new_trackedInterval_1.ndim)
    @test ([1;0] == new_trackedInterval_1.reducedDims)
    @test ([1.; 2.; 3.; 1.0] == new_trackedInterval_1.solvedVals)

    # 3d case
    Ms_2 = [[5.1;-2;2.4;.000001;.0002;.001001;;18;4;2.12;.0002;.0102;.045;;5;3.2;4.2;11;.0002;.0909;;;8.12;3.332;4.1232;11;.012002;.091239;;45.23;34;7.23;2.0;.0002;.011239;;8.127;3.2;4.3;11;.012002;.093309;;;54.1;13.332;6.1232;12.1;.12;.12239;;4.2343;5.433;94.3;8.340;.012302;.11239;;77.127;81.2;5.3;1;.02;.03], [5.836;3.122;-.031;;13.234344;42.001;.0012;;5.9;13.0022333;.0056;;;3.2;12;.0023;;56.2;1;.00233;;5.2;4.5;.013;;;5.2;4.5;.013;;3.2;12;.0023;;5.2;4.5;.013],[2;2;.002;.001;;2.1;2.1;.0021;.001;;6.23;125;.004;.0021;;;8.23453;25;.004;.0021;;2;15;.0344;.00221;;14;1.9;.004;.023021;;;8;1.9;.004;.023021;;1;1.9;.004;.023021;;4;1.9;.004;.023021]]
    errors_2 = [.0001122;.00021;.00310102]
    trackedInterval_2 = TrackedInterval([-3;4;;4;19;;-2e-16;-1e-16])
    trackedInterval_2.reducedDims = []
    trackedInterval_2.solvedVals = []
    dim_2 = 2

    expected_new_Ms_2 = [[2.70019999999999971152; 15.89019999999999832596; 0.80019999999999980034;;4.00880199999999931038; 38.00019999999999953388; 3.83900200000000069167;;48.09680000000000177351; -90.05339799999998717794; 71.84699999999999420197],[5.86699999999999999289; 13.23314399999999935176; 5.89440000000000008384;;3.19770000000000020890; 56.19767000000000223281; 5.18700000000000027711;;5.18700000000000027711; 3.19770000000000020890; 5.18700000000000027711]]
    expected_new_errors_2 = [0.00011220000000000000; 0.00021000000000000001]
    expected_new_trackedInterval_2 = [-3.;4.;;4.;19.]

    new_Ms_2, new_errors_2, new_trackedInterval_2 = reduceSolvedDim(Ms_2,errors_2,trackedInterval_2,dim_2)
    @test isapprox(expected_new_Ms_2[1],new_Ms_2[1])
    @test isapprox(expected_new_Ms_2[2],new_Ms_2[2])
    @test isapprox(expected_new_errors_2,new_errors_2)
    @test isapprox(expected_new_trackedInterval_2,new_trackedInterval_2.interval)
    @test (2 == new_trackedInterval_2.ndim)
    @test ([2] == new_trackedInterval_2.reducedDims)
    @test ([-0.00000000000000015000] == new_trackedInterval_2.solvedVals)
end