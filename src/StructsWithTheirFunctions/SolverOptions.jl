mutable struct SolverOptions
    """Settings for running interval checks, transformations, and subdivision in solvePolyRecursive.

    Parameters
    ----------
    verbose : bool
        Defaults to false. Whether or not to output progress of solving to the terminal.
    exact : bool
        Defaults to false. Whether the transformation in TransformChebInPlaceND should minimize error.
    constant_check : bool
        Defaults to true. Whether or not to run constant term check after each subdivision.
    low_dim_quadratic_check : bool
        Defaults to true. Whether or not to run quadratic check in dim 2, 3.
    all_dim_quadratic_check : bool
        Defaults to false. Whether or not to run quadratic check in dim >= 4.
    maxZoomCount : int
        Maximum number of zooms allowed before subdividing (prevents infinite infintesimal shrinking)
    level : int
        Depth of subdivision for the given interval.
    """
    verbose # = false (by default)
    exact # = false (by default)
    constant_check # = true (by default)
    low_dim_quadratic_check # = true (by default)
    all_dim_quadratic_check # = true (by default)
    maxZoomCount # = 25 (by default)
    level # = 0 (by default)
    function SolverOptions()
        #Init all the Options to default value
        new(false,false,true,true,true,25,0)
    end
end