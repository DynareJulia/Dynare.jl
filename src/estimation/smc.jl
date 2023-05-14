# Very much inspired by the Matlab implementation by Joshua Brault at 
# https://github.com/braultjosh. See also 
# 
# All references in this file unless otherwise stated relate to the textbook
# from Herbst and Schorfheide (2014) available on Frank Schorfheide's webpage

# Much of the structure of the file borrows from the SMC example used for
# Herbst and Schorfheide's (2016) textbook

using Distributions
using ProximalOperators

mutable struct Tune
    npara
    npart
    nphi
    lam
    hpdpercent
    c
    acpt
    trgt
    alp
    phi
    mu
    z
    wz
    R
    RPSD
    Rdiag
    Rchol
    Rchol2
    function Tune(nparams)
        # General
        npara = nparams                         # # of parameters
        npart = 1000                           # # of particles
        nphi  = 400                             # # of stages
        # bending coeff, lam = 1 means linear cooling schedule
        lam   = 2                               
        hpdpercent = 0.95                       # pecentage WITHIN bands
        
        # Tuning for MH algorithms
        c      = 0.5                            # initial scale cov
        acpt   = 0.25                           # initial acpt rate
        trgt   = 0.25                           # target acpt rate
        alp    = 0.9                            # Mixture weight for mixture 
                                                # proposal, 1 = only use RWMH
        
        # Create the tempering schedule
        phi = (collect(0:nphi - 1) / (nphi -1)).^lam
        mu = zeros(npara)                        # mean
        z = zeros(npart, npara)                  # centered parameters
        wz = zeros(npart, npara)                 # weighted centered parameters
        R = zeros(npara, npara)                  # covariance matrix
        RPSD = zeros(npara, npara)               # covariance matrix SPD
        Rdiag = I(npara)                         # variances
        Rchol = LowerTriangular(zeros(npara, npara)) # cholesky decomposition
        Rchol2 = I(npara)                        # standard errors
        new(npara, npart, nphi, lam, hpdpercent,
            c, acpt, trgt, alp, phi, mu, z, wz,
            R, RPSD, Rdiag, Rchol, Rchol2)
    end
end

function smc(pdraw!, logprior, loglikelihood, names, np)
    tune = Tune(np)                 

    # Matrices for storing results
    parasim = zeros(tune.nphi, tune.npart, tune.npara)     # parameter draws
    wtsim   = zeros(tune.npart, tune.nphi)                 # weights
    unwtsim   = zeros(tune.npart, tune.nphi)
    zhat    = zeros(tune.nphi)                 # normalization constant
    nresamp = 0                                # record # of iteration resampled
    csim    = zeros(tune.nphi)                 # scale parameter
    ESSsim  = zeros(tune.nphi)                 # ESS
    acptsim = zeros(tune.nphi)                 # average acceptance rate
    rsmpsim = zeros(tune.nphi)                 # 1 if re-sampled
    loglh   = zeros(tune.npart)                # log-likelihood
    logpost = zeros(tune.npart)                # log-posterior
    pdraw   = zeros(tune.npara)                # prior draws
    sig     = 0.0
    lb = - Inf
    ub = Inf

    # start by drawing from the prior distributions, with the number of draws
    # equal to the number of particles. Only accept a draw if it yields a real
    # value for the likelihood and is inside the prior bounds
    for j = 1:tune.npart
        done = false
        while !done 
            pdraw!(pdraw)
            obj = loglikelihood(pdraw)
            prior = logprior(pdraw)
            if (count(pdraw .> lb) == tune.npara 
                && count(pdraw .< ub) == tune.npara && obj > -Inf)
                parasim[1,j,:]      .= pdraw
                loglh[j]            = obj
                logpost[j]          = tune.phi[1]*obj + prior
                done                = true
            end
        end
    end
    wtsim[:, 1]    .= 1/tune.npart        # initial weight is equal weights
    zhat[1]        = sum(wtsim[:,1])     # initial normalizing constant equals 1

    # ------------------------------------------------------------------------
    # Recursion: For n=2,...,N_{\phi}
    # ------------------------------------------------------------------------
    #smctime   = tic
    #totaltime = 0
    println("SMC recursion starts ... ")

    for i = 2:tune.nphi

        #-----------------------------------
        # Step 1: Correction
        #-----------------------------------
        # incremental weights, equation (5.4)
        incwt           = exp.((tune.phi[i]-tune.phi[i-1])*loglh)
        # normalized weights (5.5), done in two steps
        unwtsim[:, i]   = wtsim[:, i-1].*incwt # numerator of normalization
        zhat[i]         = sum(unwtsim[:, i]) # denominator of normalization
        # finally normalize the weights
        wtsim[:, i]     = unwtsim[:, i]/zhat[i] 


        #-----------------------------------
        # Step 2: Selection
        #-----------------------------------
        ESS = 1/sum(wtsim[:, i].^2) # Effective sample size, (5.16)
        if (ESS < tune.npart/2) # (5.17)

            #systematic resampling
            #(id, m) = systematic_resampling(wtsim[:, i]) 
            #multinomial resampling
            id = multinomial_resampling(wtsim[:, i]) 
            @views parasim[i-1, :, :] .= parasim[i-1, id, :]
            loglh              = loglh[id]
            logpost            = logpost[id]
            # resampled weights are equal weights
            wtsim[:, i]        .= 1/tune.npart 
            nresamp            = nresamp + 1
            rsmpsim[i]         = 1 
        end

        #--------------------------------------------------------
        # (c) Mutuation
        #--------------------------------------------------------
        # Update the scaling parameter to ensure a reasonable acceptance rate.
        # See algorithm 10 on page 83.
        tune.c = tune.c*(0.95 + 0.10*exp(16*(tune.acpt-tune.trgt))/(1 +
            exp(16*(tune.acpt-tune.trgt))))
        
        # Calculate estimates of mean and variance
        @views begin
            para      = parasim[i-1, :, :]
            wght      = wtsim[:, i]
        
            # TO BE CHECKED
            mul!(tune.mu, para', wght) # mean
            tune.z .= para .- tune.mu'
            tune.wz = tune.z .* wght
            mul!(tune.R, tune.wz', tune.z)     # covariance
            prox!(tune.RPSD, IndPSD(), tune.R) # covariance PSD
        end
        tune.Rdiag   = diag(tune.RPSD) # covariance with diag elements
        tune.Rchol   .= cholesky(tune.RPSD).L  # lower triangle Cholesky 
        tune.Rchol2  = sqrt.(tune.Rdiag) 

        temp_acpt = zeros(tune.npart, 1) #initialize accpetance indicator

        # Now run a single step of RWMH-V for each particle
        for j = 1:tune.npart

            post0                   = logpost[j]
            l0                      = loglh[j]
            p0                      = para[j,:]
            # because we are at the next stage we need to account for the
            # fact that the tempering parameter is higher, update!
            post0                   = post0+(tune.phi[i]-tune.phi[i-1])*l0

            ind_acpt = 0

            # Mixture proposal if tune.alp<1, lower alp yields more draws
            # from diagonal and independence
            mix_alp = rand()
            # (a) random walk
            if mix_alp < tune.alp
                # random walk
                px     = p0 + tune.c*tune.Rchol*randn(tune.npara)
                mixsel = 1
                # (b) random walk with proposal only diagonal    
            elseif mix_alp < tune.alp + (1-tune.alp)/2
                px     = p0 + tune.c*tune.Rchol2.*randn(tune.npara)
                mixsel = 2
                # (c) independence proposal, draw starts at mean    
            else
                px     = tune.mu + tune.c*tune.Rchol*randn(tune.npara)
                mixsel = 3
            end

            # Proposal densities
            qx = logpdf_MixtureMH(px, p0, tune, mixsel)
            q0 = logpdf_MixtureMH(p0, px, tune, mixsel)

            lx = loglikelihood(px)
            prior = logprior(px)
            if any(px .< lb) || any(px .> ub)
                lx    = -1e+20
                prior = -1e+20
            end 
            postx = tune.phi[i]*lx + prior
        
            # Accept/Reject
            # this is RW, so q is canceled out
            alp = exp((postx - qx) - (post0 - q0)) 
            if rand() < alp # accept
                ind_para   = px
                ind_loglh  = lx
                ind_post   = postx
                ind_acpt   = 1
            else
                ind_para   = p0
                ind_loglh  = l0
                ind_post   = post0
                ind_acpt   = 0
            end

            # update storage matrices
            parasim[i,j,:] = ind_para
            loglh[j]       = ind_loglh
            logpost[j]     = ind_post
            temp_acpt[j,1] = ind_acpt
        end

        # update average acceptance rate
        tune.acpt = mean(temp_acpt) 
        # storage
        # scale parameter
        csim[i,:]    .= tune.c
        # ESS 
        ESSsim[i,:]  .= ESS 
        # average acceptance rate

        acptsim[i,:] .= tune.acpt 

        # The following just prints some information at the end of each stage
        if mod(i, 1) == 0
        
            para = view(parasim, i, :, :)
            wght = view(wtsim, :, i)
            mul!(tune.mu, para', wght) # mean
            tune.z .= para .- tune.mu'
            tune.wz = tune.z .* wght
            mul!(tune.R, tune.wz', tune.z)     # covariance
            prox!(tune.RPSD, IndPSD(), tune.R) # covariance PSD
            sig = sqrt.(sum(tune.z .^2 .* wght, dims=1))
            # Draws for histogram and HPD interval
            # systematic resampling
            #[id, m] = systematic_resampling(wtsim[:, i]) 
            #multinomial resampling
            id = multinomial_resampling(wght) 
            # posterior draws after resampling
            parapost = view(parasim, i, id, :)  
            # posterior probability  intervals
            # paraint = hpdint[parapost, hpdpercent, 1]  
            # time calculation
            # totaltime = totaltime + toc(smctime)
            # avgtime   = totaltime/i
            # remtime   = avgtime*(tune.nphi-i)
            # print%
            println("-----------------------------------------------")
            println(" Iteration = $i /  $(tune.nphi)")
            println("-----------------------------------------------")
            println(" phi  = $(tune.phi[i])")
            println("-----------------------------------------------")
            println("  c    = $(tune.c)")
            println("  acpt = $(tune.acpt)")
            println("  ESS  = $ESS  ($nresamp total resamples.)")
            println("-----------------------------------------------")
            #println("  time elapsed   = #5.2f\n', totaltime)
            #println("  time average   = #5.2f\n', avgtime)
            #println("  time remained  = #5.2f\n', remtime)
            #println("-----------------------------------------------")
            println("Parameter    Mean    Std     5HPD    95HPD")
            println("---------    ----    ----    ----    -----")
            for n = 1:tune.npara
                name = length(names[n]) == 2 ? "$(names[n][1]) $(names[n][2])" : names[n]
                println("$name, $(round(tune.mu[n], digits=2)), $(round(sig[n],digits=2))")#, round(paraint[1,n],digits=2), round(paraint[2,n],digits=2)])")
            end
            # smctime = tic # re-start clock
        end
    end

    # Log marginal data density approximation
    #paraint = hpdint(parapost,hpdpercent,1)  # posterior probability intervals
    println("-----------------------------------------------")
    println(" Posterior Results ")
    println("-----------------------------------------------")
    println("para      mean    std     90int")
    println("------    ----    ----     ---------")
    for n = 1:tune.npara
        name = length(names[n]) == 2 ? "$(names[n][1]) $(names[n][2])" : names[n]
        println("$name, $(round(tune.mu[n], digits=2)),$(round(sig[n], digits=2))")#, #round(paraint[1,n], digits=2), round(paraint[2,n], digits=2)])")
    end
    logmdd = sum(log(sum(unwtsim[:,2:tune.nphi], dims=1)))
    println("-----------------------------------------------")
    println(" Log marginal data density = $logmdd)")
    println("-----------------------------------------------")
    #alphaBN_loc = find(strcmp(bayestopt_.name, "alphaBN")==1)
    println("-----------------------------------------------")
    #println(" Probability of determinacy = $(sum(parapost[:,alphaBN_loc] .> 1)/tune.npart)")
    println("-----------------------------------------------")
    

    #figure_waterfall
end

"""
logpdf_MixtureMH(px, p0, tune, mixsel) 
computes the pdf of Mixture MH
    px: proposal
    p0: previous old
    tune: tuning parameters
    mixsel: mixing strategy
"""    
function logpdf_MixtureMH(px, p0, tune, mixsel)
    if mixsel == 1
        μ = p0
        Σ = tune.c^2*tune.RPSD
    elseif mixsel == 2
        μ = p0
        Σ = tune.c^2*tune.Rdiag
    elseif mixsel == 3
        μ = tune.mu
        Σ = tune.c^2*tune.RPSD
    end
    return logpdf(MvNormal(μ, Σ), px)
end    

function multinomial_resampling( w )

    np = length(w)
    cw = cumsum(w)
    
    uu = rand(np)
    indx = zeros(Int64, np)
    for i = 1:np
        
        u = uu[i]
    
        j=1
        while j <= np
            if u < cw[j]
                break
            end 
            j = j+1
        end
        indx[i] = j;
    end
    return indx
end