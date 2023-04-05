add.species.to.meta.scenario = function(newspecies, meta.scenario, taxonomy.warnings=TRUE) {
  ## Add new species to a meta.scenario. Along the way, we have to check
  ## that their classification hierarchy is consistent with the
  ## existing tree structure of the meta.scenario.
  
  taxlevels = meta.scenario$taxlevels
  Lstar = length(taxlevels)
  coding = meta.scenario$coding
  scen.classes = meta.scenario$classifications
  new.coded = matrix(NA, nrow(newspecies), Lstar)
  new.NA = rep(1, Lstar)
  for(row in seq(nrow(newspecies))) {
    for(ell in Lstar:1) {
      classification = newspecies[row, taxlevels[ell]]
      if(is.na(classification)) stop(
        sprintf(
          "Problem with entry %d (%s). Missing values are permitted in the classification only where a lower level of the classification is already known to hSSD.",
          row, newspecies[row, taxlevels[Lstar]]
        )
      )
      cint = match(classification, scen.classes[[ell]])
      if (is.na(cint)) {
        ## This classification exists and is not in the meta.scenario
        ## yet. Add it as a new classification at level ell
        scen.classes[[ell]] = c(scen.classes[[ell]], as.character(classification))
        new.coded[row, ell] = length(scen.classes[[ell]])
      } else {
        ## We have seen this classification before at level
        ## ell. Now we need to make sure that the higher levels of
        ## classification provided for the new species are
        ## consistent with ones seen before for this
        ## classification at level ell.
        
        ## Find a row with this classification in the meta.scenario and
        ## extract its coding at this and higher levels
        scenrow = match(cint, coding[,ell])
        old.higher = coding[scenrow,1:ell]
        ## Now find the coding for the higher levels provided for
        ## the new species. Don't worry about situations where the
        ## higher level is not provided.
        new.higher = sapply(
          1:ell,
          function(l) {
            cl = newspecies[row, taxlevels[l]]
            if (is.na(cl))
              NA
            else
              match(cl, scen.classes[[l]], nomatch=0)
          })
        if(!(all(old.higher[!is.na(new.higher)]==new.higher[!is.na(new.higher)]))) {
          cat("Inconsistent taxonomic classification.\n")
          cat("  Species being added is classified as:\n  ")
          print(newspecies[row,])
          cat(sprintf("  However, we already know the %s named %s. It's higher classifications are:\n  ", 
                      taxlevels[ell], newspecies[row,taxlevels[ell]]
          ))
          higher = character()
          for (i in 1:ell)
            higher = c(higher, scen.classes[[i]][old.higher[i]])
          names(higher)=taxlevels[1:ell]
          print(higher)
          stop("Exiting")
        }
        new.coded[row, 1:ell] = old.higher
        coding = rbind(coding, new.coded[row,])
        break
      }
    }
  }
  if(taxonomy.warnings) {
    print("(Re)-Classification of new species:")
    new.reclassified = sapply(
      1:Lstar,
      function(ell) scen.classes[[ell]][new.coded[,ell]]
    )
    print(new.reclassified)
  }
  meta.scenario$classifications = scen.classes
  meta.scenario$coding = coding
  meta.scenario$latin = c(meta.scenario$latin, newspecies$Latin)
  meta.scenario
}

## newdata is the log(EC50) values for tests for the new
## chemical. scenario.selection is a vector indexing the species in
## the meta-scenario in the h5file. newspecies are additional species
## to go in the scenario; they may or may not be measured. All species
## in newdata must be in either the h5 meta-scenario or newspecies or
## newspecies.testedonly.
MCMC.newchemical = function(
  newdata,
  h5file, scenario.selection, newspecies=NULL,
  newspecies.testedonly=NULL,
  continuousSSD=FALSE,
  N, burn=0, thin=1, shout=10,
  browser.on.exit = FALSE,
  detailed.output=FALSE,
  taxonomy.warnings=TRUE) {

  if(browser.on.exit)
    on.exit(browser())


  y0 = newdata$lconc
  censored = is.na(y0)
  have.censored = any(censored)
  if(have.censored) {
    n.censored = sum(censored)
    y0L = newdata$lconc.lower[censored]
    y0U = newdata$lconc.upper[censored]
    stopifnot(all(y0L<=y0U))
    stopifnot(all(is.finite(y0L)|is.finite(y0U)))
  }
  n.y0 = nrow(newdata)

  meta.scenario = h5read(h5file, "/MetaScenario")
  model = h5read(h5file, "/Model")
  q.kappa = as.vector(model$Qkappa)
  q.phi = as.vector(model$Qphi)
  N.post = h5read(h5file, "/Posterior/NoOfSamples")
  ## Could be cleverer here and read in as needed. But this will do for now.
  posterior = h5read(h5file, "/Posterior/Samples")

  ## First real job is to make the scenario by adding the new species
  ## to the relevant parts of the metascenario and then by selecting
  ## the parts of the metascenario that we want. In all the
  ## calculations which follow, we index species and taxonomic
  ## classifications with respect to the extended metascenario.
  meta.scenario = with(meta.scenario,
    list(taxlevels = NamesOfClassificationLevels,
         latin = LatinNames,
         classifications = NamesOfClassifications,
         coding = CodedClassification
         )
    )
  scenario.species = seq_along(meta.scenario$latin)[scenario.selection]
  if(!is.null(newspecies)) {
    matches = match(newspecies$Latin, meta.scenario$latin)
    if(any(is.na(matches))) {
      if(taxonomy.warnings && !all(is.na(matches))) {
          warning(sprintf("Found %d new scenario species in metascenario", sum(!is.na(matches))))
          cat("Species in newspecies which are already in the metascenario:\n")
          print(newspecies$Latin[!is.na(matches)])
      }
      meta.scenario = add.species.to.meta.scenario(
        newspecies[is.na(matches),,drop=FALSE],
        meta.scenario,
        taxonomy.warnings
        )
    } else if(taxonomy.warnings) warning("All new scenario species were already in the metascenario")
    scenario.species = c(
      scenario.species,
      match(newspecies$Latin, meta.scenario$latin)
      )
  }
  n.scenario = length(scenario.species)
  if(!continuousSSD) hc5pos = 1 + sum( ((1:n.scenario)/n.scenario) < .05 )

  if(!is.null(newspecies.testedonly)) {
    matches = match(newspecies.testedonly$Latin, meta.scenario$latin)
    if(taxonomy.warnings && !all(is.na(matches))) {
        warning(sprintf("Found %d tested-only species in metascenario", sum(!is.na(matches))))
        cat("Species in newspecies.testedonly which are already in the metascenario:\n")
        print(newspecies.testedonly$Latin[!is.na(matches)])
    }
    if(any(is.na(matches)))
       meta.scenario = add.species.to.meta.scenario(
         newspecies.testedonly[is.na(matches),,drop=FALSE],
         meta.scenario,
         taxonomy.warnings
         )
  }
  
  ## Find where measured species appear in the meta-scenario. They
  ## must do so and may or may not also be in the scenario.
  measured.species = match(newdata$latin, meta.scenario$latin)
  unique.measured.species = unique(measured.species)
  if(!all(!is.na(measured.species))) {
    cat("The following species appear in the test data but we have no taxonomic classification:\n")
    for(l in newdata$latin[is.na(measured.species)])
      cat("  ", l, "\n")
    stop("Exiting")
  }
  scenarioplus = c(
    scenario.species,
    unique.measured.species[is.na(match(unique.measured.species, scenario.species))]
    )
  index.data.in.scenarioplus = match(measured.species, scenarioplus)
  stopifnot(all(!is.na(index.data.in.scenarioplus)))

  cat(sprintf("%d species in scenario, %d species tested but not in scenario\n", n.scenario, length(scenarioplus)-n.scenario))
  
  L = model$NumberOfTaxonomicLevels
  ellstar = model$TaxonomicLevels
  meta.scenario = with(
    meta.scenario,
    list(taxlevels = taxlevels[ellstar],
         latin=latin,
         classifications = classifications[ellstar],
         coding = coding[,ellstar]
         )
    )
  stopifnot(all(sapply(model$ClassificationsUsed,
                       function(x) all(x==sort(x)))))

  ## Let's just think a little about the overall structure of the
  ## algorithm.
  ##
  ## Fundamentally, we are doing MCMC to sample the true sensitivities
  ## of all the species in the scenario to the new chemical. However,
  ## these are made from components coming from different "sources":
  ## 
  ## (i) mu and relevant beta come from the database posterior
  ## 
  ## (ii) alpha0, phi0, measured psi and measured novel beta come from
  ## the MCMC block Gibbs for the new chemical.
  ##
  ## (iii) unmeasured xi and unmeasured novel beta are sampled
  ## randomly based on sigma_xi and sigma_beta for the current
  ## selection from the database posterior.
  ##
  ## We want to make it "easy" to assemble everything from these bits.
  ## Strategy is to construct a matrix W0 which does the computation
  ## ignoring the differences between sources; this means that we have
  ## to full all the bits of beta and psi from different sources
  ## together in the right order.
  

  ## We need to know which classifications at each taxonomic level are
  ## of interest to us, i.e. which appear in the scenario. We create
  ## scenario.classifications as a list where element ell is all the
  ## values of t which appear in the scenario at level ell. In the
  ## language of the algorithms document, these are "active".
  scenario.coding = meta.scenario$coding[scenarioplus,,drop=FALSE]
  scenario.classifications = lapply(
    1:L,
    function(ell) sort(unique(scenario.coding[,ell]))
    )
  ell.scenario = rep(1:L, sapply(scenario.classifications, length))
  ## We now construct W0 with the idea that we get vector mu_0 of true
  ## sensitivities by:
  ##
  ## mu_0 = mu+alpha_0+crossprod(W0t,beta)+crossprod(W0t, psi0)
  ##
  ## where the betas (and psis) are ordered first by ell and then
  ## within ell by the meta-scenario coding (as they are in
  ## scenario.classifications).
  W0t = lapply(
    1:L,
    function(ell) as(factor(scenario.coding[,ell],
                            levels=scenario.classifications[[ell]]),
                     "sparseMatrix")
    )
  W0t = do.call("rbind", W0t)
  ## However, the betas are in an inconvenient order for what happens
  ## later. Scenario beta classifications can be: (i) relevant or
  ## novel; (ii) measured or unmeasured. We are going to reorder them
  ## as follows: RU,RM,NU,NM. To do so we need ways to select them
  ## out.

  ## Establish the coding of the measured classifications and what
  ## classifications are measured
  measured.coding = meta.scenario$coding[measured.species,,drop=FALSE]
  measured.classifications = lapply(
    1:L,
    function(ell) sort(unique(measured.coding[,ell]))
    )

  ## We can now classify the active (ell, t) in two ways. Firstly, are
  ## they "measured" or not?
  scenario.classifications.measured = lapply(
    1:L,
    function(ell) !is.na(match(scenario.classifications[[ell]],
                               measured.classifications[[ell]]))
    )
  ## Turn into one long vector of TRUE/FALSE
  scM = unlist(scenario.classifications.measured)
  n.psi.M = sum(scM)
  n.psi.U = sum(!scM)
  ell.psi.M = ell.scenario[scM]
  ell.psi.U = ell.scenario[!scM]
  ## Secondly, does beta_{\ell t} fail to appear in the posterior from the
  ## database? This is "novelty" in the language of the algorithms
  ## document.
  scenario.classifications.novel = lapply(
    1:L,
    function(ell) is.na(match(scenario.classifications[[ell]],
                              model$ClassificationsUsed[[ell]]))
    )
  ## Turn into one long vector of TRUE/FALSE
  scN = unlist(scenario.classifications.novel)
  ## Make the four categories of scenario classification
  scRU = (!scN) & (!scM)
  scRM = (!scN) & scM
  scNU = scN & (!scM)
  scNM = scN & scM
  n.beta.NU = sum(scNU)
  ell.beta.NM = ell.scenario[scNM]
  ell.beta.NU = ell.scenario[scNU]
  ## Order the columns of W0beta as desired. Within each of the four
  ## groups they are still in the original order.
  W0betat = rbind(
    W0t[scRU,,drop=FALSE],
    W0t[scRM,,drop=FALSE],
    W0t[scNU,,drop=FALSE],
    W0t[scNM,,drop=FALSE]
    )
  ## Psi is more simple: just unmeasured (first) or measured (second)
  W0psit = rbind(
    W0t[!scM,,drop=FALSE],
    W0t[scM,,drop=FALSE]
    )
  
  ## We also want to know which of the beta variables in the data
  ## posterior are "relevant":
  model.classifications.relevant = lapply(
    1:L,
    function(ell) !is.na(match(model$ClassificationsUsed[[ell]],
                               scenario.classifications[[ell]]))
    )

  ## Extract the relevant betas from the posterior and bind together
  ## to a single matrix
  post.beta = lapply(
    1:L,
    function(ell)
      posterior$beta[[ell]][,model.classifications.relevant[[ell]],drop=FALSE]
  )
  stopifnot(all(
    sapply(post.beta, ncol)+sapply(scenario.classifications.novel, sum) ==
    sapply(scenario.classifications, length))
    )
  post.beta = do.call("cbind", post.beta)

  stopifnot(sum(scRU)+sum(scRM)==ncol(post.beta))
  post.beta.U = post.beta[,scRU[scRU|scRM],drop=FALSE]
  post.beta.M = post.beta[,scRM[scRU|scRM],drop=FALSE]


  ## Now we compute the matrices W_1 and W_2 for use in the new
  ## chemical MCMC. This bit all relates only to measured (ell, t)
  Wt.betapsi = lapply(
    1:L,
    function(ell) as(factor(measured.coding[,ell],
                            levels=measured.classifications[[ell]]),
                     "sparseMatrix")
  )

  measured.classifications.novel = lapply(
    1:L,
    function(ell) is.na(match(measured.classifications[[ell]],
                              model$ClassificationsUsed[[ell]]))
    )
  Wt.beta.novel = lapply(
    1:L,
    function(ell) Wt.betapsi[[ell]][measured.classifications.novel[[ell]],,drop=FALSE]
    )
  stopifnot(all(ell.beta.NM==rep(1:L, sapply(Wt.beta.novel, nrow))))
  Wt.beta.novel = do.call("rbind", Wt.beta.novel)
  stopifnot(nrow(Wt.beta.novel)==sum(scNM))


  Wt.beta.relevant = lapply(
    1:L,
    function(ell) Wt.betapsi[[ell]][!measured.classifications.novel[[ell]],,drop=FALSE]
    )
  Wt.beta.relevant = do.call("rbind", Wt.beta.relevant)
  stopifnot(nrow(Wt.beta.relevant)==sum(scRM))
  
  W1t = rbind(
    Matrix(1, 1, n.y0),  # mu
    Wt.beta.relevant   # measured beta in the database
    )

  ell.psi =  rep(1:L, sapply(Wt.betapsi, nrow))
  Wt.betapsi = do.call("rbind", Wt.betapsi)
  W2t = rbind(
    Matrix(1, 1, n.y0),  # alpha_0
    Wt.beta.novel,      # measured beta not in the database
    Wt.betapsi          # measured psi
    )
  W2t = Matrix(W2t, sparse=TRUE)
  n.vartheta = nrow(W2t)
  beta.lookup = 1+seq_along(ell.beta.NM)
  psi.lookup = 1+length(ell.beta.NM)+seq_along(ell.psi)



  
  ## Initial values
  if(have.censored) 
    y0[censored] = ifelse(
        is.finite(y0L),
        ifelse(is.finite(y0U), (y0L+y0U)/2, y0L),
        y0U
        )
  kappa0 = rep(1, n.y0)
  phi0 = sd(y0)
  post.row = sample(N.post, 1)
  beta.RU = post.beta.U[post.row,]
  beta.RM = post.beta.M[post.row,]
  mu = posterior$mu[post.row]
  sigma.epsilon = posterior$sigma$epsilon[post.row]
  sigma.alpha = posterior$sigma$alpha[post.row]
  sigmas.beta = posterior$sigma$beta[post.row,]
  sigmas.xi = posterior$sigma$xi[post.row,]
  nu.kappa = posterior$nu$kappa[post.row]
  nu.phi = posterior$nu$phi[post.row]

  ## Storage for results
  hc5.out = numeric(N)
  accept.out = logical(N)
  if (detailed.output)
    mu0.out = matrix(0, N, length(scenarioplus))
  
  cat("Starting main loop\n")
  
  for(t in seq(-burn+1,N)) {
    if(!is.null(shout))
      if (t%%shout==0)
        cat(sprintf("Iteration %d\n", t))
    anyaccepted = FALSE
    for(dummy in 1:thin) {
      ## Do the mixed model parameter update
      sigmatilde = c(
        sigma.alpha,
        sigmas.beta[ell.beta.NM],
        sigmas.xi[ell.psi]*phi0
        )

      rtilde = sigma.epsilon/sqrt(kappa0)
      varthetastar = rnorm(n.vartheta, 0, sigmatilde)
      estar = rnorm(n.y0, 0, rtilde)
      SighWtRnegh = Diagonal(x=sigmatilde) %*% W2t %*% Diagonal(x=1/rtilde)
      y = y0 - crossprod(W1t, c(mu, beta.RM))
      RHS = SighWtRnegh %*% ((y-crossprod(W2t, varthetastar)-estar)/rtilde)
      Cstarchol = Cholesky(tcrossprod(SighWtRnegh), Imult=1)
      varthetatilde = sigmatilde*solve(Cstarchol, RHS)
      vartheta = varthetatilde+varthetastar

      alpha0 = vartheta[1]
      beta.NM = vartheta[beta.lookup]
      psi.M = vartheta[psi.lookup]

      ## Sample the kappa values for the new chemical data
      epsilon0 = as(y - crossprod(W2t, vartheta), "numeric")
      kappa0 = rgamma(
        n.y0,
        (nu.kappa+1)/2,
        (nu.kappa-q.kappa+(epsilon0/sigma.epsilon)^2)/2
        )

      ## Sample phi_0
      phi0 = 1/sqrt(rgamma(
        1,
        (nu.phi+n.psi.M)/2,
        (nu.phi-q.phi+sum(psi.M^2/sigmas.xi[ell.psi]^2))/2
        ))

      ## Sample the novel unmeasured beta and psi
      beta.NU = rnorm(n.beta.NU, 0, sigmas.beta[ell.beta.NU])
      psi.U = rnorm(n.psi.U, 0, sigmas.xi[ell.psi.U]*phi0)

      ## Compute the true sensitivities to the new chemical of all the
      ## species in the scenario
      beta = c(beta.RU, beta.RM, beta.NU, beta.NM)
      psi = c(psi.U, psi.M)

      mu0 = mu+alpha0+crossprod(W0betat, beta)+crossprod(W0psit, psi)

      ## Sample values for the censored new chemical data (if any) 
      if(have.censored) 
        y0[censored] = rcensnorm(
            n.censored,
            y0L,
            y0U,
            mu0[index.data.in.scenarioplus[censored]],
            sigma.epsilon/sqrt(kappa0[censored])
            )
      
      ## Compute the HC5
      if (continuousSSD)
        hc5 = quantile(mu0[1:n.scenario], .05)
      else
        hc5 = sort(mu0[1:n.scenario])[hc5pos]

      ## Propose a new sample from the posterior from the database
      post.row = sample(N.post, 1)
      ## Compute the acceptance ratio
      prop.beta.RU = post.beta.U[post.row,]
      prop.beta.RM = post.beta.M[post.row,]
      prop.mu = posterior$mu[post.row]
      prop.sigma.epsilon = posterior$sigma$epsilon[post.row]
      prop.sigma.alpha = posterior$sigma$alpha[post.row]
      prop.sigmas.beta = posterior$sigma$beta[post.row,]
      prop.sigmas.xi = posterior$sigma$xi[post.row,]
      prop.nu.kappa = posterior$nu$kappa[post.row]
      prop.nu.phi = posterior$nu$phi[post.row]

      prop.beta = c(prop.beta.RU, prop.beta.RM, beta.NU, beta.NM)
      prop.mu0 = prop.mu+alpha0+crossprod(W0betat, prop.beta)+crossprod(W0psit, psi)
      
      lr1 = sum(
        dnorm(y0,
              prop.mu0[index.data.in.scenarioplus],
              prop.sigma.epsilon/sqrt(kappa0),
              log=TRUE) -
        dnorm(y0,
              mu0[index.data.in.scenarioplus],
              sigma.epsilon/sqrt(kappa0),
              log=TRUE)
        )
      lr2 = dnorm(alpha0, 0, prop.sigma.alpha, log=TRUE)-
        dnorm(alpha0, 0, sigma.alpha, log=TRUE)
      lr3 = dgamma(1/phi0^2, prop.nu.phi/2, (prop.nu.phi-q.phi)/2, log=TRUE)-
        dgamma(1/phi0^2, nu.phi/2, (nu.phi-q.phi)/2, log=TRUE)
      lr4 = sum(
        dgamma(kappa0, prop.nu.kappa/2, (prop.nu.kappa-q.kappa)/2, log=TRUE)-
        dgamma(kappa0, nu.kappa/2, (nu.kappa-q.kappa)/2, log=TRUE)
        )
      lr5 = sum(
        dnorm(psi.M, 0, phi0*prop.sigmas.xi[ell.psi.M], log=TRUE)-
        dnorm(psi.M, 0, phi0*sigmas.xi[ell.psi.M], log=TRUE)
        )
      lr6 = sum(
        dnorm(beta.NM, 0, prop.sigmas.beta[ell.beta.NM], log=TRUE)-
        dnorm(beta.NM, 0, sigmas.beta[ell.beta.NM], log=TRUE)
        )
      lrsum = lr1+lr2+lr3+lr4+lr5+lr6

      ## Accept or reject
      accept = log(runif(1))<=lrsum

      if(accept) {
        beta.RU = prop.beta.RU
        beta.RM = prop.beta.RM
        mu = prop.mu
        sigma.epsilon = prop.sigma.epsilon
        sigma.alpha = prop.sigma.alpha
        sigmas.beta = prop.sigmas.beta
        sigmas.xi = prop.sigmas.xi
        nu.kappa = prop.nu.kappa
        nu.phi = prop.nu.phi
      }

      anyaccepted = anyaccepted | accept

    }
    if(t<0) next
    ## Save results
    hc5.out[t] = hc5
    accept.out[t] = anyaccepted
    if(detailed.output) 
      mu0.out[t,] = as.numeric(mu0)


  }
  ret = list(hc5=hc5.out, accept=accept.out, p=(hc5pos-1)/n.scenario)
  if(detailed.output) {
    ret$n.scenario = n.scenario
    colnames(mu0.out) = meta.scenario$latin[scenarioplus]
    ret$mu0  = mu0.out
  }
  
  ret

}
