include("../../Julia-Rootfinding/yroots/ChebyshevApproximator.jl")

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