function quadraticCheck2D(test_coeff,tol)
    if ndims(test_coeff) != 2
        return false
    end
    interval = [[-1, -1], [1, 1]]
# PRINT STATEMENTS TO TEST C
    c = 1.0*fill(0,6)
    shape = reverse(size(test_coeff))
    c[1] = test_coeff[1,1]
    if shape[1] > 1
        c[2] = test_coeff[1,2]
    end
    if shape[2] > 1
        c[3] = test_coeff[2,1]
    end
    if shape[1] > 2
        c[4] = test_coeff[1,3]
    end
    if shape[1] > 1 && shape[2] > 1
        c[5] = test_coeff[2,2]
    end
    if shape[2] > 2
        c[6] = test_coeff[3,1]
    end

    other_sum = sum(abs.(test_coeff)) - sum(abs.(c)) + tol
    k0 = c[1]-c[4]-c[6]
    k3 = 2*c[4]
    k5 = 2*c[6]
    function eval_func(x,y)
        return k0 + (c[2] + k3 * x + c[5] * y) * x  + (c[3] + k5 * y) * y
    end

    det = 16 * c[4] * c[6] - c[5]^2
    if det != 0
        int_x = (c[3] * c[5] - 4 * c[2] * c[6]) / det
        int_y = (c[2] * c[5] - 4 * c[3] * c[4]) / det
    else                      # det is zero,
        int_x = Inf
        int_y = Inf
    end

    min_satisfied, max_satisfied = false,false
    #Check all the corners
    eval = eval_func(interval[1][1], interval[1][2])
    min_satisfied = min_satisfied || eval < other_sum
    max_satisfied = max_satisfied || eval > -other_sum
    if min_satisfied && max_satisfied
        return false
    end

    eval = eval_func(interval[2][1], interval[1][2])
    min_satisfied = min_satisfied || eval < other_sum
    max_satisfied = max_satisfied || eval > -other_sum
    if min_satisfied && max_satisfied
        return false
    end

    eval = eval_func(interval[1][1], interval[2][2])
    min_satisfied = min_satisfied || eval < other_sum
    max_satisfied = max_satisfied || eval > -other_sum
    if min_satisfied && max_satisfied
        return false
    end

    eval = eval_func(interval[2][1], interval[2][2])
    min_satisfied = min_satisfied || eval < other_sum
    max_satisfied = max_satisfied || eval > -other_sum
    if min_satisfied && max_satisfied
        return false
    end

    #Check the x constant boundaries
    #The partial with respect to y is zero
    #Dy:  c4x + 4c5y = -c2 =>   y = (-c2-c4x)/(4c5)
    if c[6] != 0
        cc5 = 4 * c[6]
        x = interval[1][1]
        y = -(c[3] + c[5]*x)/cc5
        if interval[1][2] < y < interval[2][2]
            eval = eval_func(x,y)
            min_satisfied = min_satisfied || eval < other_sum
            max_satisfied = max_satisfied || eval > -other_sum
            if min_satisfied && max_satisfied
                return false
            end
        end
        x = interval[2][1]
        y = -(c[3] + c[5]*x)/cc5
        if interval[1][2] < y < interval[2][2]
            eval = eval_func(x,y)
            min_satisfied = min_satisfied || eval < other_sum
            max_satisfied = max_satisfied || eval > -other_sum
            if min_satisfied && max_satisfied
                return false
            end
        end
    end

    #Check the y constant boundaries
    #The partial with respect to x is zero
    #Dx: 4c3x +  c4y = -c1  =>  x = (-c1-c4y)/(4c3)
    if c[4] != 0
        cc3 = 4*c[4]
        y = interval[1][2]
        x = -(c[2] + c[5]*y)/cc3
        if interval[1][1] < x < interval[2][1]
            eval = eval_func(x,y)
            min_satisfied = min_satisfied || eval < other_sum
            max_satisfied = max_satisfied || eval > -other_sum
            if min_satisfied && max_satisfied
                return false
            end
        end

        y = interval[2][2]
        x = -(c[2] + c[5]*y)/cc3
        if interval[1][1] < x < interval[2][1]
            eval = eval_func(x,y)
            min_satisfied = min_satisfied || eval < other_sum
            max_satisfied = max_satisfied || eval > -other_sum
            if min_satisfied && max_satisfied
                return false
            end
        end
    end

    #Check the interior value
    if interval[1][1] < int_x < interval[2][1] && interval[1][2] < int_y < interval[2][2]
        eval = eval_func(int_x,int_y)
        min_satisfied = min_satisfied || eval < other_sum
        max_satisfied = max_satisfied || eval > -other_sum
        if min_satisfied && max_satisfied
            return false
        end
    end

    return true
end