using Combinatorics
using AxisArrays: AxisArray
using DualNumbers
using Distributions
using FiniteDiff
using LinearAlgebra
import NLsolve
using NonlinearSolve
using Roots
using Tasmanian
using Printf

export DDSG_test1, DDSG_test2

const DDSGInx    = Tuple{Int,Vararg{Int}}
const DDSGInxMap = Dict{DDSGInx, Int}

# Define the DDSG struct with updated field names and added attributes X0 and Y0
mutable struct DDSG
    dim::Int                           # Dimension of the grid
    dof::Int                           # Degrees of freedom for grid output
    l_min::Int                         # Minimum sparse grid level
    l_max::Int                         # Maximum sparse grid level
    k_max::Int                         # Maximum combination length   
    order::Int                         # Polynomial order for local basis
    rule::String                       # Rule for constructing the sparse grid
    domain::Matrix{Float64}            # Domain as a matrix
    centroid::Vector{Float64}          # Center of the domain as a vector
    lookup::Vector{DDSGInxMap}         # Vector of dictionaries to hold grid structures
    rlookup::Vector{DDSGInx}           # Vector of tuples for reverse lookup
    grid::Vector{Any}                  # Vector to hold grid structures
    coeff::Vector{Float64}             # Vector to hold grid DD coefficients
    grid_points::Vector{Int}           # Vector to hold grid points
    X0::Vector{Float64}                # Anchor point X0
    Y0::Vector{Float64}                # Values at anchor point
    coeff0::Float64                    # Coefficient for the zeroth DD order component function
    is_ddsg::Bool                      # is ddsg grid

    # Define constructor for the struct
    function DDSG(; 
        dim::Int, 
        dof::Int, 
        l_min::Int, 
        l_max::Int, 
        k_max::Int, 
        order::Int = 1,
        rule::String = "localp",
        domain::Union{Matrix{Float64},Nothing}=nothing, 
        lookup::Vector{DDSGInxMap} = Vector{DDSGInxMap}(), 
        rlookup::Vector{DDSGInx} = Vector{DDSGInx}(),  
        grid::Vector{Any} = Vector{Any}(), 
        coeff::Vector{Float64} = Vector{Float64}(), 
        grid_points::Vector{Int} = Vector{Int}(), 
        X0::Vector{Float64} = Vector{Float64}(), 
        Y0::Vector{Float64} = Vector{Float64}(), 
        coeff0::Float64 = 0.0,
        is_ddsg::Bool = true
    )

        if isnothing(domain)
            domain = zeros(Float64,dim,2)
            domain[:,2].=1.0
        end

        # Set center as the midpoint of each dimension in the domain
        centroid = vec((domain[:, 1] + domain[:, 2]) / 2.0)

        # Resize lookup to length of k_max
        resize!(lookup, k_max)

        # Initialize empty dictionaries for the resized lookup
        for i in 1:k_max
            lookup[i] = DDSGInxMap()
        end

        if dim==k_max
            is_ddsg = false
        end

        return new(
            dim, 
            dof, 
            l_min, 
            l_max, 
            k_max, 
            order, rule, domain, centroid, 
            lookup,rlookup, grid, coeff, grid_points, 
            X0, Y0, coeff0, is_ddsg)
    end

    # Function to create a copy of the DDSG instance
    function DDSG(other::DDSG)

        # Deep copy all attributes
        dim_copy         = deepcopy(other.dim)
        dof_copy         = deepcopy(other.dof)
        l_min_copy       = deepcopy(other.l_min)
        l_max_copy       = deepcopy(other.l_max)
        order_copy       = deepcopy(other.order)
        rule_copy        = deepcopy(other.rule)
        k_max_copy       = deepcopy(other.k_max)    
        domain_copy      = deepcopy(other.domain)
        centroid_copy    = deepcopy(other.centroid)
        coeff_copy       = deepcopy(other.coeff)
        grid_points_copy = deepcopy(other.grid_points)
        lookup_copy      = deepcopy(other.lookup)
        rlookup_copy     = deepcopy(other.rlookup)
        X0_copy          = deepcopy(other.X0)
        Y0_copy          = deepcopy(other.Y0)
        coeff0_copy      = deepcopy(other.coeff0)
        is_ddsg_copy     = deepcopy(other.is_ddsg)     

        # Specific deep copy for grid using copyGrid function
        grid_copy = Vector{Any}(undef, length(other.grid))
        for i in 1:length(other.grid)
            grid_copy[i] = copyGrid(other.grid[i])  # Use copyGrid to deep copy each grid
        end

        # Create a new DDSG instance with the copied values
        return new(
            dim_copy,
            dof_copy, 
            l_min_copy, 
            l_max_copy, 
            k_max_copy,    
            order_copy, rule_copy, domain_copy, centroid_copy,
            lookup_copy, rlookup_copy,grid_copy, coeff_copy, grid_points_copy, 
            X0_copy, Y0_copy, coeff0_copy,is_ddsg_copy)
    end
end

# Constructor method to initialize DDSG and set up grid and center
function DDSG_init!(self::DDSG)

    if !self.is_ddsg

        k = self.dim
        u = Tuple(1:k)
        domain_loc = self.domain
                
        # Initialize a Tasmanian sparse grid for the current combination
        grid_u = Tasmanian.TasmanianSG(k, self.dof, self.l_min)
        makeLocalPolynomialGrid!(grid_u, order=self.order, rule=self.rule)
        setDomainTransform!(grid_u, domain_loc)  # Match grid to domain dimensions
        #enableAcceleration(grid_u,"none")

        #print(grid_u)
        #print(domain_loc)
        
        # Assign default values to vectors and map vector indices using lookup
        push!(self.grid, grid_u)
        push!(self.coeff, 1.0)
        push!(self.grid_points, 0.0)
        push!(self.rlookup, u)
        self.lookup[k][u] = length(self.grid)
    else
        # Loop over each possible combination length up to k_max
        for k in 1:self.k_max
            for u in combinations(1:self.dim, k)
                # Slice the domain by the selected combination and convert to a tuple for dictionary key
                domain_loc = self.domain[u, :]
                u = Tuple(u)
                
                # Initialize a Tasmanian sparse grid for the current combination
                grid_u = Tasmanian.TasmanianSG(k, self.dof, self.l_min)
                makeLocalPolynomialGrid!(grid_u, order=self.order, rule=self.rule)
                setDomainTransform!(grid_u, domain_loc)  # Match grid to domain dimensions
                
                # Assign default values to vectors and map vector indices using lookup
                push!(self.grid, grid_u)
                push!(self.coeff, 0.0)
                push!(self.grid_points, 0.0)
                push!(self.rlookup, u)
                self.lookup[k][u] = length(self.grid)
            end
        end

    end
end

# Method to print DDSG instance details
function DDSG_print(self::DDSG)
    println("\n=== DDSG Object Details START ===")
    println("Dimension (dim): ", self.dim)
    println("Degrees of freedom (dof): ", self.dof)
    println("Level min (l_min): ", self.l_min)
    println("Level max (l_max): ", self.l_max)
    println("Polynomial order (order): ", self.order)
    println("Rule (rule): ", self.rule)
    println("Domain (domain): ", self.domain)
    println("Domain centroid (centroid): ", self.centroid)
    println("Maximum expansion order (k_max): ", self.k_max)
    println("Anchor point (X0): ", self.X0)
    println("Function values at anchor point (Y0): ", self.Y0)
    println("DD coefficient for 0-order (coeff0): ", self.coeff0)
    println("Is DDSG (is_ddsg): ", self.is_ddsg)
    println("Grid details:")
    for k in 1:length(self.lookup)
        println("- Order (k): ", k)
        for (u, inx) in self.lookup[k]
            println("-- Index (inx): ", inx)
            println("-- Index (u): ", u)
            println("-- DD coefficient (coeff[lookup[k][u]]): ", self.coeff[inx])
            println("-- Grid points (grid_points[lookup[k][u]]): ", self.grid_points[inx])
            println(self.grid[inx])
        end
    end
    println("=== DDSG Object Details END ===")
end


# Function to print matrix with fixed number of decimal places
function ppnice(matrix, digits=4)
    for row in eachrow(matrix)
        for element in row
            @printf("%.*f  ", digits, element)
        end
        println()  # Newline after each row
    end
end


# Function to build the DDSG instance, iterating over grids and updating points
function DDSG_build!(
    self::DDSG; 
    F::Function, 
    X0::Vector{Float64},
    refinement_tol::Float64,
    refinement_type::String="classic",
    scale_corr_vec::Union{Vector{Float64},Nothing}=nothing)
    
    #t_build = 0

    self.X0 = X0
    self.Y0 = vec(F(reshape(X0, :, 1)))

    #println("\n=== DDSG Build Start ===")
    for k in 1:length(self.lookup)
       # println("-- Order (k): $k")

        for (u, inx) in self.lookup[k]

            u_inx  = collect(u)
            println("    Index (u): $u")

            # Load the grid
            grid = self.grid[inx]
            self.grid_points[inx] = 0

            # Refinement based on surplus coefficients
            for l in (self.l_min):self.l_max
     
                if getNumNeeded(grid) < 1
                    break
                end

                #println("l $l l_max $(self.l_max) points-$(getNumNeeded(grid)) ")

                # Get and prepare needed points for the grid
                X = getNeededPoints(grid)

                # println(" x_points ",X)

                N = size(X,2)

                X_full = repeat(self.X0, 1, N)

                X_full[u_inx,:] = X

                # for row in eachrow(X_full)
                #     println(join(round.(row, digits=3), "  "))
                # end

                # readline()

                #println(" -> to be evaluated X $(round.(X_full,digits=3))")

                #println("-- Level $l  size $(size(X_full))  normX $(norm(X_full)) normY $(norm(Y_val)) ")

                #println("    Points needed refinement level $l: $N")
                #t_build = t_build - time()
                Y_val = F(X_full)
                #t_build = t_build + time()


                loadNeededPoints!(grid, Y_val)

                #println("-+ Level $l  size $(size(X_full))  meanX $(mean(X_full)) meanY $(mean(Y_val)) ")

                self.grid_points[inx]+=N

                #println(grid)
                if isnothing(scale_corr_vec)
                    # output=-1 refine a selects which output to use for refinement sequence and local polynomial grids accept -1 to indicate all outputs
                    setSurplusRefinement!(grid, refinement_tol, output=-1, refinement_type=refinement_type)
                else 
                    scale_corr_mat = repeat(scale_corr_vec, 1, getNumLoaded(grid))
                    # if size(scale_corr_mat, 2) < getNumLoaded(grid)
                    #     scale_corr_mat = repeat(Float64.(scale_corr), 1, getNumLoaded(grid))
                    # end
                    setSurplusRefinement!(grid, refinement_tol, output=-1, refinement_type=refinement_type, scale_correction=scale_corr_mat')
                end
             
            end
        end

        #println(" t_build. inside ... $t_build")
    end

    if self.is_ddsg 

        for k in 1:length(self.lookup)
            for (u, _) in self.lookup[k]
                len_u = length(u)
                #println("- u=$u")
                
                for i in len_u:-1:1
                    V = combinations(u,i)
                    len_u_less_len_v = len_u - i
                    for v in V
                        inx = self.lookup[i][Tuple(v)]
                        #println("  v=$v")
                        self.coeff[inx]+= (-1)^len_u_less_len_v
                    end
                end

                # for v=()
                len_u_less_len_v = len_u - 0
                self.coeff0+= (-1)^len_u_less_len_v

            end
        end

        # for u=()
        len_u_less_len_v = 0 - 0
        self.coeff0+= (-1)^len_u_less_len_v

    else
        self.coeff0=0
        self.coeff[end]=1
    end

end

function DDSG_evaluate!(ddsg::DDSG; Y::AbstractVecOrMat{Float64}, X::AbstractVecOrMat{Float64})

  

    is_batch = false
    if X isa AbstractMatrix
        is_batch = true
    end

    if !ddsg.is_ddsg

        if is_batch
            evaluateBatch!(Y, ddsg.grid[end], X) 
            #Y. = evaluateBatch( ddsg.grid[end], X) 
            #println("----evaluateBatch! eval ",norm(X),"  ",norm(Y))
            #readline()
        else
            Y .= evaluate(ddsg.grid[end],vec(X))
        end

    else
        if is_batch
            for j in 1:size(Y, 2)
                Y[:,j] .= ddsg.Y0
            end
        else
            Y .= ddsg.Y0
        end

        Y .= Y .* ddsg.coeff0

        for (inx, u) in enumerate(ddsg.rlookup)

            # Load the grid
            grid = ddsg.grid[inx]
            u_inx = collect(u)

            if is_batch
                X_partial = X[u_inx,:]
                Y_buffer  = zeros(size(Y))
                evaluateBatch!(Y_buffer, grid, X_partial) 
            else
                X_partial = X[u_inx]
                Y_buffer  = evaluate(grid,vec(X_partial))
            end

            Y .= Y .+ Y_buffer * ddsg.coeff[inx]
        end

    end

    # println(X," ",size(X))
    # println(Y," ",size(Y))
    # println(ddsg.domain," ",size(ddsg.domain))


    # Y_copy = copy(Y)  # Create a copy of the original matrix

    # # a = minimum(Y[end,:])
    # # b = minimum(X[end,:])
    # # print(a," ", b, "\n")
   
    # for i in 1:(size(Y, 1) - 1)  # Iterate over all rows except the last one
    #     Y[i, :] .= clamp.(Y[i, :], ddsg.domain[i, 1], ddsg.domain[i, 2])  # Apply clipping row-wise
    # end



        
    # if any(Y .!= Y_copy) 

    #     println("Y_base: \n",Y_copy,size(Y_copy))
    #     println("Y: \n",Y,size(Y))

    #     readline()
    # end



end

function DDSG_evaluate(ddsg::DDSG; X::AbstractVecOrMat{Float64})
    Y = zeros(Float64,ddsg.dof,size(X,2))
    DDSG_evaluate!(ddsg,Y=Y,X=X)
    return Y
end

function DDSG_differentiate!(ddsg::DDSG; J::Matrix{Float64}, x::Vector{Float64})

    if !ddsg.is_ddsg
        differentiate!(J,ddsg.grid[end],x)
    else
        for (inx, u) in enumerate(ddsg.rlookup)
            grid = ddsg.grid[inx]
            u_inx = collect(u)
            J_buffer = similar(J, length(u_inx), ddsg.dof)
            x_partial = x[u_inx]

            differentiate!(J_buffer,grid,x_partial)
            J[u_inx,:] .+= J_buffer * ddsg.coeff[inx]
        end
    end
end

function DDSG_nodes(ddsg::DDSG)

    N = sum(ddsg.grid_points)

    # Initialize an empty matrix to store the cumulative horizontal concatenation
    X_loc_total = Float64[]  # Start with an empty matrix

    for (inx, u) in enumerate(ddsg.rlookup)
        # Load the grid
        grid = ddsg.grid[inx]
        u_inx = collect(u)
        X_partial = getPoints(grid)

        # Repeat ddsg.X0 to match the size of X_partial
        X_loc = repeat(ddsg.X0, 1, size(X_partial, 2))
        X_loc[u_inx, :] .= X_partial

        # Concatenate X_loc horizontally with X_loc_total
        if length(X_loc_total) == 0
            X_loc_total = X_loc
        else
            X_loc_total = hcat(X_loc_total, X_loc)
        end
    end

    return X_loc_total  # Return the concatenated result
end

function  DDSG_values(ddsg::DDSG)

    N = sum(ddsg.grid_points)

    # Initialize an empty matrix to store the cumulative horizontal concatenation
    X_loc_total = Float64[]  # Start with an empty matrix

    for (inx, u) in enumerate(ddsg.rlookup)
        # Load the grid
        grid = ddsg.grid[inx]
        u_inx = collect(u)
        X_partial = getLoadedValues(grid)

        X_loc_total = hcat(X_loc_total, X_loc)
    end

    return X_loc_total  # Return the concatenated result
end

function DDSG_test_unit(;dim,k_max, c,
    dof=10, l_min=2, l_max=8, tol=1e-6, num_points=1000)

    F_poly = (X_in; c, dof) -> begin
        if ndims(X_in) == 1
            X_in = reshape(X_in, :, 1)
        end
        result = Matrix{Float64}(undef, dof, size(X_in,2))
        for j in 1:size(X_in,2)
            for i in 1:dof 
                temp = sin.(  X_in[:, j] )
                result[i,j] =  sum(temp)^Int(c)

            end
        end
        return result
    end


    F = X->F_poly(X,c=c,dof=dof)

    if true
        grid = DDSG(dim=dim, dof=dof,l_min=l_min,l_max=l_max, k_max=k_max)
        DDSG_init!(grid)

        DDSG_build!(grid,F=F,X0=grid.centroid,refinement_tol=tol) 
        Random.seed!(1234)
        X_sample = Matrix{Float64}(undef, dim, num_points)
        for i in 1:num_points
            X_sample[:, i] .= grid.domain[:, 1] .+ (grid.domain[:, 2] - grid.domain[:, 1]) .* rand(dim)
        end

        Y = DDSG_evaluate(grid,X=X_sample)
        num_points = sum(grid.grid_points)
    else
        domain = zeros(dim,2)
        domain[:,2].=1.0
        grid = TasmanianSG(dim, dof, l_min)
        makeLocalPolynomialGrid!(grid, order=1, rule="localp")
        setDomainTransform!(grid, domain)  # Match grid to domain dimensions

        # Refinement based on surplus coefficients
        for l in (l_min):l_max
            if getNumNeeded(grid) < 1
                break
            end
            X = getNeededPoints(grid)
            N = size(X,2)
            Y_val = F(X)
            loadNeededPoints!(grid, Y_val)
            println("-+ Level $l  size $(size(X))  meanX $(mean(X)) meanY $(mean(Y_val)) ")
            setSurplusRefinement!(grid, tol, refinement_type="classic")
        end
        Random.seed!(1234)
        X_sample = Matrix{Float64}(undef, dim, num_points)
        for i in 1:num_points
            X_sample[:, i] .= domain[:, 1] .+ (domain[:, 2] - domain[:, 1]) .* rand(dim)
        end

        Y = Tasmanian.evaluateBatch(grid,X_sample)
        num_points = getNumLoaded(grid)
    end

    Y_exact = F(X_sample)

    return norm(Y - Y_exact)/norm(Y_exact), num_points
end

function DDSG_test1(;dim,k_max,l_vec,c)

    err_vec = []
    pts_vec = []
    tim_vec = []

    for l_max in l_vec

        if k_max == -1
            k_max_loc = dim
        else
            k_max_loc = k_max
        end

        time_taken = @elapsed begin
            (err,pts) = DDSG_test_unit(
                dim=dim, 
                k_max=k_max_loc, 
                c=c, 
                dof=10,
                l_min=l_max, # Non adaptive
                l_max=l_max)
        end

        push!(tim_vec,time_taken)
        push!(err_vec,err)
        push!(pts_vec,pts)
    end

    println("Dim=$dim k_max=$k_max c=$c" )
    println("tim_vec=$tim_vec")
    println("err_vec=$err_vec")
    println("pts_vec=$pts_vec")
    println("l_vec=$l_vec")

    println("")
end


function DDSG_test2(;l,k_max,dim_vec,c)

    err_vec = []
    pts_vec = []
    tim_vec = []

    for dim in dim_vec

        if k_max == -1
            k_max_loc = dim
        else
            k_max_loc = k_max
        end

        time_taken = @elapsed begin
            (err,pts) = DDSG_test_unit(
                dim=dim, 
                k_max=k_max_loc, 
                c=c, 
                dof=10,
                l_min=l, # Non adaptive
                l_max=l)
        end

        push!(tim_vec,time_taken)
        push!(err_vec,err)
        push!(pts_vec,pts)
    end

    println("l=$l k_max=$k_max c=$c" )
    println("tim_vec=$tim_vec")
    println("err_vec=$err_vec")
    println("pts_vec=$pts_vec")  
    println("dim_vec=$dim_vec")

    println("")
end

function ex2()
    dim =  2       
    outs = 1       
    iDepth = 2
    tol = 1e-6   
    K = 7  # max refinement steps
    tsg = Tasmanian.TasmanianSG(dim,outs,iDepth)
    which_basis = 1 #1= linear basis functions -> Check the manual for other options
    Tasmanian.makeLocalPolynomialGrid!(tsg,order=which_basis,rule="localp")

    negbox = x-> x*2 - 1
    
    # sparse grid points from that object
    spPoints = getPoints(tsg)
    
    # test fun 
    tfun(x,y) = exp(-x^2) * cos(y)

    Random.seed!(2)
    N = 1000
    randPnts = negbox.(rand(2,N))
    # truth
    truth = [tfun(randPnts[1,i], randPnts[2,i]) for i in 1:N]

    # values on sparse grid
    spVals = [tfun(spPoints[1,i], spPoints[2,i]) for i in 1:size(spPoints,2)]
    # load points needed for such values
    loadNeededPoints!(tsg,spVals)

    # evaluate interpolation
    res = evaluateBatch(tsg,randPnts)

    numpoints = size(spPoints,2)

    @info("error on initial grid:    $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

    # refinefment loop
    #        anim = @animate for k in 1:K
    for k in 1:K
        setSurplusRefinement!(tsg, tol, refinement_type="classic")
        if Tasmanian.getNumNeeded(tsg) > 0
            spPoints = getNeededPoints(tsg)   # additional set of points required after refinement
            spVals = [tfun(spPoints[1,i], spPoints[2,i])::Float64 for i in 1:size(spPoints,2)]
            # load points needed for such values
            loadNeededPoints!(tsg, spVals)
            numpoints =+ size(spPoints,2)
            
            # evaluate interpolation
            res = evaluateBatch(tsg,randPnts)
        pred = Tasmanian.evaluateBatch(tsg,Array(spPoints))  
            @info("refinement level $k error: $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")
            
            # plot
            zerone = (-1.1,1.1)
            #=
        plot(scatter(spPoints[1, :],spPoints[2, :],title="level $k grid:\n $numpoints points",m=(:black,1,:+),aspect_ratio=:equal,xlims=zerone,ylims=zerone),
             scatter3d(randPnts[1, :],randPnts[2, :],res[:],title="max error: $(round(maximum(abs,res .- truth),digits = 5))",m=(:red,1),xlims=zerone,ylims=zerone,zlims=(0,1)),
             scatter3d(spPoints[1, :],spPoints[2, :],pred[:],title="grid prediction",m=(:red,1,0.2),zlims=(0,1),xlims=zerone,ylims=zerone,zgrid=:black),
             layout=(1,3),leg=false
        )
            =#
        end
    end
    #if save 
    #    gif(anim,joinpath(dirname(@__FILE__),"ex2.gif"),fps=1)
    #end
    return nothing
end
