################################################################################
#                           Time iteration step                                #
################################################################################


function ti_step(grid, pol_guess, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, y)

    # Get the points that require function values
    aPoints1 = Tasmanian.getNeededPoints(grid)
    # Get the number of points that require function values
    aNumAdd = Tasmanian.getNumNeeded(grid)

    # Array for intermediate update step
    polInt = zeros(length(system_variables), aNumAdd)
    let state
        fnl(x) = sysOfEqs(x, y, state, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables)
        
        # Time Iteration step
        for ii1 in 1:1#aNumAdd
            @views state = aPoints1[:, ii1]
            @views pol = pol_guess[:, ii1]
            res = Dynare.nlsolve(fnl, pol, method = :robust_trust_region, show_trace = false)
            polInt[:, ii1] .= res.zero
            @show polInt[:, ii1]
        end
#        ln(-1)
    end

    # Add the new function values to grid1
    Tasmanian.loadNeededPoints!(grid, polInt)

    return grid
end


################################################################################
#                             Grid refinement                                  #
################################################################################

function refine(grid, nPols, scaleCorr, surplThreshold, dimRef, typeRefinement)

    # Get the points that require function values
    aNumLoad = Tasmanian.getNumLoaded(grid)
    # Scaling to only allow for those policies that are supposed to
    # determine the refinement process (in this case only the capital policies)
    scaleCorrMat = zeros(nPols, aNumLoad )
    scaleCorrMat .= scaleCorr
    
    # Refine the grid based on the surplus coefficients
    Tasmanian.setSurplusRefinement!(grid, surplThreshold, iOutput=dimRef, sCriteria=typeRefinement, llfScaleCorrection=scaleCorrMat)
    if Tasmanian.getNumNeeded(grid) > 0

	# Get the new points and the number of points
	nwpts = Matrix(Tasmanian.getNeededPoints(grid))
	aNumNew = Tasmanian.getNumNeeded(grid)
        
	# We assign (for now) function values through interpolation#
	pol_guess = zeros(aNumNew, nPols)
	pol_guess = Tasmanian.evaluateBatch(grid, nwpts)
    else
        
	pol_guess = []

    end
    return grid, pol_guess
end


################################################################################
#                        New grid construction                                 #
################################################################################

function fresh_grid(gridDim, gridOut, gridDepth, gridDomain, gridOrder, gridRule)

    # Generate the grid structure
    grid = Tasmanian.TasmanianSG(gridDim, gridOut, gridDepth)
    Tasmanian.makeLocalPolynomialGrid!(grid, iOrder = gridOrder, sRule = gridRule)
    # Transform the domain
    Tasmanian.setDomainTransform!(grid, gridDomain)

    return grid
end


################################################################################
#                Checking convergence and updating policies                    #
################################################################################

function policy_update(gridOld, gridNew, nPols)

    # Get the points and the number of points from grid1
    aPoints2 = Tasmanian.getPoints(gridNew)
    aNumTot = Tasmanian.getNumPoints(gridNew)

    # Evaluate the grid points on both grid structures
    polGuessTr1 = Tasmanian.evaluateBatch(gridNew, Matrix(aPoints2))
    polGuessTr0 = Tasmanian.evaluateBatch(gridOld, Matrix(aPoints2))

    # 1) Compute the Sup-Norm

    metricAux = zeros(nPols)

    for imet in 1:nPols
        @views metricAux[imet] = maximum(abs.(polGuessTr0[:, imet] - polGuessTr1[:, imet]))
    end
    metricSup = maximum(metricAux)

    # 2) Compute the L2-Norm

    metricL2 = 0.0

    for imetL2 in 1:nPols
        @views metricL2 += sum(abs.(polGuessTr0[:, imetL2] - polGuessTr1[:, imetL2]).^2)
    end
    
    metricL2 = (metricL2/(aNumTot*nPols))^0.5

    metric = min(metricL2, metricSup)

    # Now update pol_guess and grid

    polGuess = zeros(nPols, aNumTot)
    @show size(polGuess)
    @show size(polGuessTr1)
    for iupd in 1:nPols
        @views polGuess[:,iupd] = 0.5*polGuessTr0[:,iupd] + 0.5*polGuessTr1[:,iupd]
    end

    ### !!! to be implemented !!!
    gridOld = Tasmanian.copyGrid(gridNew)

    return metric, polGuess, gridOld
end

################################################################################
#                               Grid storage                                   #
################################################################################

function save_grid(grid, iter)

    grid.write(data_location_smooth + "grid_iter_$i.txt")

    return
end
