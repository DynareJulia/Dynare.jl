# All references in this file unless otherwise stated relate to the textbook from Herbst and Schorfheide (2014) available on Frank
# Schorfheide's webpage

# Much of the structure of the file borrows from the SMC example used for
# Herbst and Schorfheide's (2016) textbook
#addpath('./Functions');
#addpath('./dynare/5.2/matlab');


#dynare Smets_Wouters_2007_45
#[dataset_, dataset_info, ~, ~, M_, options_, oo_, estim_params_,bayestopt_, bounds] = dynare_estimation_init({}, M_.dname, 0, M_, options_, oo_, estim_params_, bayestopt_);
#fcn = @evaluate_likelihood2;

function prior_draw!(pdraw, ep)
    for (i, p) in enumerate(ep.prior)
        pdraw[i] = rand(p, 1)
    end
end

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
    function Tune(pshape)
        # General
        npara = length(bayestopt_.pshape)       # # of parameters
        npart = 10000                           # # of particles
        nphi  = 400                             # # of stages
        lam   = 2                               # # bending coeff, lam = 1 means linear cooling schedule
        hpdpercent = 0.95                       # pecentage WITHIN bands
        
        # Tuning for MH algorithms
        c      = 0.5                            # initial scale cov
        acpt   = 0.25                           # initial acpt rate
        trgt   = 0.25                           # target acpt rate
        alp    = 0.9                            # Mixture weight for mixture proposal, 1 = only use RWMH
        
        # Create the tempering schedule
        phi = (collect(0:nphi - 1) / (nphi -1)).^lam
        new(npara, npart, nphi, lam, hdpercent,
            c, acpt, trgt, alp, phi)
    end
end

function smc(;context = context)
    tune = Tune(pshape)                     

    # Matrices for storing results
    parasim = zeros(tune.nphi, tune.npart, tune.npara)     # parameter draws
    wtsim   = zeros(tune.npart, tune.nphi)                 # weights
    zhat    = zeros(tune.nphi)                             # normalization constant
    nresamp = 0                                            # record # of iteration resampled
    csim    = zeros(tune.nphi)                             # scale parameter
    ESSsim  = zeros(tune.nphi)                             # ESS
    acptsim = zeros(tune.nphi)                             # average acceptance rate
    rsmpsim = zeros(tune.nphi)                             # 1 if re-sampled
    loglh   = zeros(tune.npart)                            # log-likelihood
    logpost = zeros(tune.npart)                            # log-posterior
    pdraw   = zeros(tune.npara)                            # prior draws

    # start by drawing from the prior distributions, with the number of draws
    # equal to the number of particles. Only accept a draw if it yields a real
    # value for the likelihood and is inside the prior bounds
    for j = 1:tune.npart
        done = false
        while done 
            prior_draw!(pdraw)
            (obj, prior)            = feval(fcn, pdraw, M_, estim_params_, oo_, options_, bayestopt_)
            if count(pdraw .> bounds.lb) == tune.npara && count(pdraw < bounds.ub) == tune.npara && obj > -Inf
                parasim[1,j,:]      .= pdraw
                loglh[j]            = obj
                logpost[j]          = tune.phi[1]*obj + prior
                done                = true
            end
        end
    end

    wtsim[:, 1]    = 1/tune.npart                          # initial weight is equal weights
    zhat[1]        = sum(wtsim(:,1))                       # initial normalizing constant equals 1

    # ------------------------------------------------------------------------
    # Recursion: For n=2,...,N_{\phi}
    # ------------------------------------------------------------------------
    #smctime   = tic
    #totaltime = 0
    @show "SMC recursion starts ... "

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
        ESS = 1/sum(wtsim(:, i).^2) # Effective sample size, (5.16)
        if (ESS < tune.npart/2) # (5.17)
            #(id, m) = systematic_resampling(wtsim[:, i]) #systematic resampling
            (id, m) = multinomial_resampling(wtsim[:, i]) #multinomial resampling
            @views parasim[i-1, :, :] .= parasim[i-1, id, :]
            loglh              = loglh[id]
            logpost            = logpost[id]
            wtsim[:, i]        = 1/tune.npart # resampled weights are equal weights
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
        para      = squeeze(parasim[i-1, :, :])
        wght      = repeat(wtsim[:, i], 1, tune.npara)
        
        tune.mu      = sum(para.*wght) # mean
        z            = (para - repeat(tune.mu, tune.npart, 1))
        #tune.R       = (z.*wght)'*z       # covariance
        tune.R       = nearestSPD((z.*wght)'*z)
        tune.Rdiag   = diag(diag(tune.R)) # covariance with diag elements
        tune.Rchol   = chol(tune.R)
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
        mix_alp = rand
        # (a) random walk
        if mix_alp < tune.alp
            # random walk
            px     = p0 + tune.c*tune.Rchol*randn(tune.npara)
            mixsel = 1
            # (b) random walk with proposal only diagonal    
        elseif mix_alp < tune.alp + (1-tune.alp)/2
            px     = p0 + tune.c*tune.Rchol2*randn(tune.npara)
            mixsel = 2
            # (c) independence proposal, draw starts at mean    
        else
            px     = tune.mu + tune.c*tune.Rchol*randn(tune.npara)
            mixsel = 3
        end

        # Proposal densities
        qx = logpdf_MixtureMH(px, p0, tune, mixsel)
        q0 = logpdf_MixtureMH(p0, px, tune, mixsel)

        
        (lx, prior)           = feval(fcn, px, M_, estim_params_, oo_, options_, bayestopt_)
        if any(px < bounds.lb) || any(px > bounds.ub)
            lx    = -1e+20
            prior = -1e+20
        end
        postx                   = tune.phi(i)*lx + prior
        
        # Accept/Reject
        alp = exp((postx - qx) - (post0 - q0)) # this is RW, so q is canceled out
        if rand < alp # accept
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

    tune.acpt = mean(temp_acpt) # update average acceptance rate
    # storage
    csim[i,:]    = tune.c # scale parameter
    ESSsim[i,:]  = ESS # ESS
    acptsim[i,:] = tune.acpt # average acceptance rate


    # The following just prints some information at the end of each stage
    if mod(i, 1) == 0
        
        para = squeeze(parasim[i, :, :])
        wght = repmat(wtsim[:, i], 1, tune.npara)
        mu  = sum(para.*wght)
        sig = sum((para - repeat(mu, tune.npart)).^2 .* wght)
        sig = (sqrt(sig))
        # Draws for histogram and HPD interval
        #[id, m] = systematic_resampling(wtsim[:, i]) # systematic resampling
        (id, m) = multinomial_resampling(wtsim[:, i]) #multinomial resampling
        parapost = squeeze(parasim[i, id, :])  # posterior draws after resampling
        paraint = hpdint[parapost, hpdpercent, 1]  # posterior probability intervals
        # time calculation
        # totaltime = totaltime + toc(smctime)
        # avgtime   = totaltime/i
        # remtime   = avgtime*(tune.nphi-i)
        # print
        @show "-----------------------------------------------"
        @show " Iteration = $i /  $(tune.nphi)"
        @show "-----------------------------------------------"
        @show " phi  = $(tune.phi[i])"
        @show "-----------------------------------------------"
        @show "  c    = $(tune.c)"
        @show "  acpt = $(tune.acpt)"
        @show "  ESS  = $ESS  ($nresamp total resamples.)"
        @show "-----------------------------------------------"
        #@show "  time elapsed   = #5.2f\n', totaltime)
        #@show "  time average   = #5.2f\n', avgtime)
        #@show "  time remained  = #5.2f\n', remtime)
        #@show "-----------------------------------------------"
        @show "Parameter    Mean    Std     5HPD    95HPD"
        @show "---------    ----    ----    ----    -----"
        for n = 1:tune.npara
            @show "$([bayestopt_.name[n], round(mu[n],2),round(sig[n],2), round(paraint[1,n],2), round(paraint[2,n],2)]"
        end
        # smctime = tic # re-start clock
    end

    # Log marginal data density approximation
    paraint = hpdint(parapost,hpdpercent,1)  # posterior probability intervals
    @show "-----------------------------------------------"
    @show " Posterior Results "
    @show "-----------------------------------------------"
    @show "para      mean    std     90int"
    @show "------    ----    ----     ---------"
    for n = 1:tune.npara
        @show "$([bayestopt_.name[n], round(mu[n],2),round(sig[n],2), round(paraint[1,n],2), round(paraint[2,n],2)]"
    end
    logmdd = sum(log(sum(unwtsim(:,2:tune.nphi))))
    @show "-----------------------------------------------"
    @show " Log marginal data density = $logmdd)"
    @show "-----------------------------------------------"
    alphaBN_loc = find(strcmp(bayestopt_.name, "alphaBN")==1)
    @show "-----------------------------------------------"
    @show " Probability of determinacy = $(sum(parapost[:,alphaBN_loc] .> 1)/tune.npart)"
    @show "-----------------------------------------------"
    

    #figure_waterfall
end

