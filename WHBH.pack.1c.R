# LT 20/1/2018
#
# 20/01/2018
# First release candidate of between host - within host simulation functions
# RC0
#
# 12/02/2018
# fixed the scheduler so no event can happen after tmax
# 
# Tajima's D function for Seedy simulations modified
# meansnps, the nucleotide diversity function of seedy is broken
# corrected the calculations of pi in WHBH by a small factor of (n-1)/n to account 
# for the diagonal of the distance matrix
# implemented an alignment function for whbh data objects absed on seedy 2dna

require(deSolve)
require(igraph)

#' @title Stochastic simulation of evolutionary dynamics during SIR epidemics 
#'
#' @description
#' \code{sim.epi} simulates an SIR disease outbreak with explicit representation
#' of within host dynamics...
#' @param nH A positive integer. Number of hosts in the population at the
#' beginning of the outbreak. One host is infected, the remainder are
#' susceptible. 
#' @param tmax A numerical scalar. Simulation will stop when \code{tmax}
#' is reached if the outbreak has not finished yet. 
#' @param sample.n A positive integer. Number of individuals to sample around
#' peak infection. 
#' @param seed A integer to seed the random number generator.
#' @param verbose A logical scalar. Should information on the unfolding of the
#' outbreak be printed during the simulation ? 
#' @param mu a positive numerical scalar. The natural death rate of hosts in
#' the absence of infection.
#' @param bcst a positive numerical scalar. Per capita transmissibility rate
#' of the pathogen.
#' @param vcst a positive numerical scalar. 
#' @param mut a positive numerical scalar. The per generation mutation rate.
#' @param mut.cutoff 
#' @param r 
#' @param P0 
#' @param rM 
#' @param M0 
#' @param delta 
#'
#' @return a whbhsim object, containing the results of the simulation:
#' $mutation.tree. A vector of positive integers representing the mutations
#' ...
#' 
#' @export
#'
#' @examples
sim.epi = function(nH=100, tmax=1e3, sample.n = nH/10, seed = NULL, verbose=FALSE,
                   mu = 1e-3, bcst=0.3, vcst=0.01,
                   mut = 3e-3, mut.cutoff=1e2, r=1, P0=1e-9, rM = 0.5, M0 = 1e-7, delta = 1.5) {
  # nH is the number of hosts, all susceptibles at the begining except one infected
  # tmax: maximum simulation time
  # sample.n: number of hosts sampled around peak infection time
  ## between host parameters
  # mu: natural death rate
  # bcst: relative transmissibility
  # vcst: relative virulence
  # mut.cutoff: relative prevalence cutoff under which mutations are ignored in within host dynamics
  # r: instrinsic growth rate of pathogen
  # P0: ratio of pathogen's initial population size over maximum pop. size
  # rM: intrinsic growth rate of the immune system
  # M0: immune system's ratio of initial population size over maximum pop. size
  # delta: per capita attack rate of immune system
  # DELTA SHOULD BE BIGGER than r
  
  start.t = proc.time()
  if (!is.null(seed))
    set.seed(seed)
  
  ## work out within host dynamics
  
  # time at peak infection 
  # according to exp-logistic approx.
  tp = -1/rM*log(M0/(1-M0)*(delta/r-1))
  
  # time to end of infection
  # exp-logistic approx
  te = tp - r*tp/(r-delta)
  
  ## ODE simulation
  f = function(t,state,par) {
    with(as.list(c(state, par)), {
      dP = r * (1-P) * P - delta * M * P
      dM = rM * (1-M) * M
      return(list(c(dP,dM)))
    })
  }
  
  par = list(r=r, rM=rM, delta=delta)
  
  Ps = ode(c(P=P0, M=M0), seq(0,2*te, length.out = 1000), f, par)
  
  ## work out simulated values for tp, te, etc.
  Pp = max(Ps[,"P"])
  tp = Ps[which.max(Ps[,"P"]), "time"]
  te = Ps[min(which(Ps[,"P"]<P0)),"time"]
  
  # redo ode simulation until te
  Ps = ode(c(P=P0, M=M0), seq(0,te, length.out = 1000), f, par)
  
  # tabulate the integral
  intP = c(0,cumsum(Ps[,"P"])[-nrow(Ps)]*(Ps[2,"time"]- Ps[1,"time"]))
  
  # create interpolation of integral of patogen poplation size
  intP.fun = approxfun(Ps[,"time"], intP, rule=2)
  
  ################################################################################
  ## states
  ################################################################################
  
  # number of susceptibles in the population
  S = nH-1
  
  # number of infected
  I = 1
  
  #list of hosts
  #hosts = c()
  
  # first infected host
  #hosts[[1]] = host()
  
  # vector of infection time
  I.t = 0
  
  # recovered
  R=0
  
  states = c(t=0, S=S, I=I, R=R)
  states.sav = states
  
  
  # pathogen's strains tracking
  patho = 0
  names(patho)=0
  
  # association between hosts and strains
  I.patho = 1
  
  # mutation info
  I.muts = c()
  
  # transmission tree
  trans.tree = list()
  trans.tree[[1]] = c(0, 1, 0, NA, NA)
  
  # association between infected host and transmission tree
  I.trans = 1
  
  ################################################################################
  ##
  
  ################################################################################
  #
  ## rates function
  #
  ################################################################################
  
  
  # this function builds a rate vector, as a function of time and states
  rates = function(t, t0, states, params, sum=TRUE) {
    # we consider 4 types of events: transmission, death of host per virulence,
    # death of host per natural death, recovery
    # becuase the process in inhomogenous, we return the integral of the rate
    # with time
    # t: current time
    # t0: time at which last event occured.
    
    I = states["I"]
    S = states["S"]
    R = states["R"]
    
    iPs = intP.fun(t-I.t) - intP.fun(t0-I.t)
    transmit = bcst * iPs * S / (S+I+R)
    
    # mortality by virulence
    virul = vcst * iPs 
    
    # natural mortality, constant
    mort = mu * (t-t0) 
    
    # recovery rate, constant
    # time to extinction exp/logis
    #recov = 1/te * (t - t0)
    # count recovery time from time of infection
    recov = 1/te * (t - pmax(t0, I.t))
    
    rates = rbind(transmit, virul, mort, recov)
    
    if (sum) return(sum(rates))
    rates  
  }
  
  cleanInfected = function(ihost) {
    I.t <<- I.t[-ihost]
    I.patho <<- I.patho[-ihost]
    I.muts <<- I.muts[-ihost]
    I.trans <<- I.trans[-ihost]
  }
  
  execute.event = function(ihost, ievent, t, states) {
    
    # save time of event
    states["t"] = t
    
    if (ievent == 0) {
      # transmission
      if (verbose) cat("host ", ihost, " just transmitted a mysterious disease to a new host\n")
      states["S"] = states["S"] - 1
      states["I"] = states["I"] + 1
      I.t[states["I"]] <<- t
      
      # get transmitted strain, update new host strain
      I.patho[states["I"]] <<- sample.strain(ihost)
      
      # record transmission in transmission tree
      trans.tree[[length(trans.tree) + 1]] <<- c(I.trans[ihost], I.patho[states["I"]], t, 0 , 0)
      I.trans[states["I"]] <<- length(trans.tree)
      
    } else if (ievent == 1) {
      # killed
      if (verbose) cat ("host ", ihost, " was killed by a mysterious virus\n")
      states["I"] = states["I"] - 1
      
      #record time and cause of death
      trans.tree[[I.trans[ihost]]][4:5] <<- c(t,1)
      
      cleanInfected(ihost)
      
    } else if (ievent == 2) {
      # dead
      if (verbose) cat ("host ", ihost, " died from natural death\n")
      states["I"] = states["I"] - 1
      
      #record time and cause of death
      trans.tree[[I.trans[ihost]]][4:5] <<- c(t,2)
      
      cleanInfected(ihost)
      
    } else if (ievent == 3) {
      # recovery
      if (verbose) cat("host ", ihost, " recovered from a long disease\n")
      states["I"] = states["I"] - 1
      
      #record time and cause of death
      trans.tree[[I.trans[ihost]]][4:5] <<- c(t,3)
      
      cleanInfected(ihost)
      states["R"] = states["R"] + 1
    }
    states
  }
  
  ## mutation business
  
  # generate mutations 
  getMutations = function(mut.cutoff) {
    t = 0
    imut=2
    muts = c('0'=1)
    tree = 0
    
    # 1e4 is the arbitrary limit. means this strain had 1e-4 probability of being transmitted
    while(t<tp && muts[length(muts)] < mut.cutoff) {
      p = runif(1)
      t.mut = uniroot(function(s) r*mut/P0*(intP.fun(s) - intP.fun(t)) + log(p), c(t,tmax), tol=1e-8)$root
      
      # work out who the father is
      ps.cum = cumsum(c(0,1/muts))
      tree[imut] = findInterval(runif(1, max=max(ps.cum)), ps.cum)
      
      # keep population size (exponential approx., could do better)
      muts[imut] = exp(r*t.mut)
      names(muts)[imut] = t.mut
      imut = imut+1
      t = t.mut  
    }
    rbind(muts, tree)
  }
  
  sample.strain = function(ihost) {
    # sample one strain from the many strains present in the host
    # build list of host strains if necessary and stores in I.muts
    # return the index of sampled strain  
    
    # do we already have mutation information for this host ?
    if (length(I.muts)< ihost || is.null(I.muts[[ihost]])) {
      if (verbose) cat("Bulding genetic information for host", ihost,"\n")
      I.muts[[ihost]] <<- getMutations(mut.cutoff)
    }
    
    # mutation stuff: is the transmitted pathogen a mutant ?
    muts = I.muts[[ihost]]
    # draw strain/lineage
    # relative abundances of strains
    ps=cumsum(c(0, 1/muts[1,1:findInterval(t-I.t[ihost], as.numeric(colnames(muts)))]))
    ps = ps/ps[length(ps)]
    
    # draw random strain
    p = runif(1)
    strain = findInterval(p, ps)
    
    # if selected strain id the wild, we are done. return the patho node
    if (strain == 1) return(I.patho[[ihost]])
    
    ## if a mutant is transmitted, rebuild the tree, update patho
    
    # go up the tree until the root is reached
    lineage=strain
    
    while (strain > 1) {
      # parent of current strain
      strain = muts[2,strain]
      lineage = c(strain, lineage)
    }
    
    # remove the wild (root) from lineage
    lineage = lineage[-1]
    
    if (length(lineage)>1) {
      if (verbose) cat("a secondary mutant has been transmited\n")
      
      # new pathogens, set appropriate lineage
      patho[length(patho)+1:length(lineage)] <<- c(I.patho[ihost], length(patho)+1:(length(lineage)-1))
      # update time-stamps
      names(patho)[length(patho)-((length(lineage)-1):0)] <<- as.numeric(colnames(muts)[lineage]) + I.t[ihost]
      
    } else {
      # new pathogen, with the wild as father
      patho[length(patho)+1] <<- I.patho[ihost]
      names(patho)[length(patho)] <<- as.numeric(colnames(muts)[lineage]) + I.t[ihost]
    }
    
    # return the added strain
    length(patho)
    
  }
  
  # sampling  by scared epidemiologists
  sample.epidemic = function(states, n=nH/10) {
    
    if (verbose) cat("\n\n SAMPLING of epidemic by WHO\n\n")  
    # choose n individuals among infected
    obs = sample(states["I"], min(n, states["I"]))
    
    # for each sample individual get sampled genome
    # problem here is that we could sample a strain which is not in patho
    # because it has never been transmitted
    # first order approximation: if a non wild is picked up, assume it has not
    # been transmitted and update patho accordingly
    strains = sapply(obs, sample.strain)
    
    # return smapled strains
    strains
  }
  
  ################################################################################
  #
  ## main loop
  #
  ################################################################################
  
  t = 0
  sampled = FALSE
  sampled.2 = FALSE
  sampled.time = 0
  sampled.strain=c()
  sampled.strain.2=c()
  
  while (t < tmax) {
    
    ## do we still have infected ppl ?
    if (states["I"] <= 0) break;
    
    ## draw time of next event
    # random probability
    p = runif(1)
    
    # check next event occurs before we overshoot tmax:
    if( rates(tmax,t, states, params) + log(p) < 0 ) break;
    
    t.next = uniroot(function(s) rates(t+s, t, states, params) + log(p), c(0,tmax-t), tol=1e-6)
    if (t.next$root == 0) {
      # FIXME test whether uniroot gives a sensical answer
      t.next = uniroot(function(s) rates(t+s, t, states, params) + log(p), c(0,tmax - t), tol=1e-10)
    } 
    
    t.new = t + t.next$root
    
    ## event treatment
    # find what event occured
    event = which(rmultinom(1,1,prob=rates(t.new, t, states, params, sum=FALSE))==1) - 1
    
    t = t.new
    
    # host index
    ihost = event %/% 4 + 1
    
    # type of event
    ievent = event %% 4
    
    states = execute.event(ihost, ievent, t, states)
    
    # check for smapling of disease
    # sample at peak infection
    # thtas when S/S+I+R = 1/R0
    if (!sampled && states["S"]/(sum(states)) < 1/intP.fun(te)/bcst) {
      sampled = TRUE
      sampled.t = t
      sampled.strain = sample.epidemic(states, n=sample.n)
      sampled.I = states["I"]
    }
    
    # double sampling, when I is half what it was at first sampling
    if (sampled && !sampled.2 && states["I"] < sampled.I/5) {
      sampled.2 = TRUE
      sampled.strain.2 = sample.epidemic(states, n=sample.n)
    }
    
    # save result
    states.sav = rbind(states.sav, states)
  }
  if (verbose) {cat("\n"); print(proc.time() - start.t); cat("\n\n")}
  
  row.names(states.sav) = 1:nrow(states.sav)
  whbhsim = list(mutation.tree = patho,
             trans.tree = matrix(unlist(trans.tree), byrow=T, ncol=5,
                                 dimnames=list(1:length(trans.tree), c("tree", "patho", "t0", "te", "cause"))),
             sampled.strain = sampled.strain,
             sampled.time = ifelse(sampled, sampled.t, -1),
             sampled.strain.2 = sampled.strain.2,
             bh.dynamics = states.sav)
  attr(whbhsim, "class") = "whbhsim"
  
  whbhsim
}

plot.whbhsim = function (sim, ...) {
  if (!is.null(sim$bh.dynamics))
    matplot(sim$bh.dynamics[,1], sim$bh.dynamics[,-1], type="l")
}


# plot tree
plottree = function(tree, subset=NA, cols=NA) {
  require(igraph)
  
  # ensure root has no dad
  if (length(tree) > 1)
    gtree=graph_from_data_frame(data.frame(from=2:length(tree),to=tree[-1]), vertices=data.frame(vertex=1:length(tree)))
  else
    return()
  
  if (length(cols)==1) {
    cols = c(1,rep(2,length(tree)-1))
    if (length(subset) > 1)
      cols[subset]=3
  }
  #  plot(gtree, layout=layout_as_tree(gtree, circular=T), vertex.size=3, vertex.label=NA,vertex.color=cols, edge.arrow.size=.3, edge.width=2)
  plot(gtree, vertex.size=3, vertex.label=NA,vertex.color=cols, edge.arrow.size=.3, edge.width=2)
  
  
}

Dn = function(whbhsim) {
  
  # function to calculate distance between nodes of pathogen tree
  tree.dist = function(i, j) {
    if (i < j) return(tree.dist(j, i))
    if (i == j) return(0)
    if (i == 1) return(root.dist(j))
    tree.dist(whbhsim$mutation.tree[i],j) + 1
  }
  
  # function to calculate distance from node to root
  root.dist = function(i) {
    if (i==0) return(-1)
    root.dist(whbhsim$mutation.tree[i]) + 1
  }
  
  # nucleotide diversity
  D = 0
  
  # number of seggregating sites (Watterson estimator)
  W = 0
  
  strains = whbhsim$sampled.strain
  
  # check if got sampled data 
  if (is.null(strains) || length(unique(strains)) <= 1) {
    res=c(res, pi=NA, W=NA, D=NA)
    next
  }

  n = length(strains)
  
  # iterate on alleles
  for (i in 1:n) {
    dists = sapply(i:n, function(j) tree.dist(strains[i], strains[j]))
    D = D + sum(dists)
    W = max(W, max(dists))
  }

  return(c(D=D * 2 / n / (n-1), W = W/sum(1/(1:(n-1)))))
  
}

# lots of Dn:

getDns = function(n=100) {
  replicate(n, {cat("New sim\n"); repeat { sim = sim.epi(10000, sample.n=1000, mu=0, mut=1e-3); if(sim$sampled.time > 0) break; }; Dn(sim) })
}

TajimaD = function(whbhsim) {
  
  # function to calculate distance between nodes of pathogen tree
  tree.dist = function(i, j) {
    if (i < j) return(tree.dist(j, i))
    if (i == j) return(0)
    #if (i == 1) return(root.dist(j))
    tree.dist(whbhsim$mutation.tree[i],j) + 1
  }

  # this function returns the path between two nodes as
  # a vector of complex numbers. Each complex number represent one edge
  # in the format node1 + i node2, with node1 > node2
  tree.dist.2 = function(i, j) {
    if (i < j) return(tree.dist.2(j, i))
    if (i == j) return(c())
    dad.i = whbhsim$mutation.tree[i]
    c(complex(1,i, dad.i), tree.dist.2(dad.i,j))
  }
  
  # function to calculate distance from node to root
  # useless, commented
  #root.dist = function(i) {
  #  if (i==0) return(-1)
  #  root.dist(whbhsim$mutation.tree[i]) + 1
  #}
  
  res= c()
  
  #dirty as hell, hack to check how peak infection and end of infection compare
  for (strains in list(whbhsim$sampled.strain, whbhsim$sampled.strain.2)) {
    

  # nucleotide diversity
  p = 0
  
  # list of edges (= subtree) needed to join every strain to each other
  # number of seggregating sites
  Ss = c()
  
  # check if got sampled data 
  if (is.null(strains) || length(unique(strains)) <= 1) {
    res = c(res,pi=NA,W=NA,D=NA)
    next
  }
  
  N = length(strains)
  
  # get frequency of alleles
  alleles = table(strains)
  fs = alleles/sum(alleles)
  n = length(fs)
  strains = as.numeric(names(fs))
  
  # iterate on alleles
  for (i in 1:(n-1)) {
    edges = lapply((i+1):n, function(j) tree.dist.2(strains[i], strains[j]))
    dists = sapply(edges, length)
    Ss = unique(c(Ss, unlist(edges)))
    p = p + fs[i] * sum(dists * fs[(i+1):n])
  }
  
  # number of segregating sites
  S = length(Ss)
  
  # average nucleotide diversity
  # the factor 2 comes from the fact that I iterated only on the upper triangle
  # of the distance matrix
  # the correction factor comes from the absence of the diagonal in the
  # calculations of Tajima (1983)
  pi = 2 * N / (N-1) *unname(p)
  
  # sampling corrections and standard deviation of difference
  # see https://en.wikipedia.org/wiki/Tajima%27s_D or
  # Tajima, F. (1989) Statistical method for testing the neutral mutation
  # hypothesis by DNA polymorphism. Genetics, 123, 595--595.
  
  a1 = sum(1/(1:(N-1)))
  a2 = sum(1/(1:(N-1))^2)
  b1 = (N+1)/3/(N-1)
  b2 = 2 * (N^2 + N + 3)/9/N/(N-1)
  c1 = b1 - 1/a1
  c2 = b2 - (N+2)/a1/N + a2/a1^2
  e1 = c1 / a1
  e2 = c2 / (a1^2 + a2)

  # Watterson estimator
  W = S/a1 
  
  # Tajima's D
  D = (pi - W)/sqrt(e1 * S + e2 * S * (S-1))
  
  res = c(res, c(pi = pi, W = W, D = D))

  }
  
  return(res)
  
}

# simulation study on behavior of Tajima's D:

simulate = function(N=1000) {

  rs = c(1, 1.5, 2)
  deltas = c(1.1, 1.5, 2)
  nHs = c(1000, 1e4)
  
  res = array(NA, dim = c(N, 6, length(rs), length(deltas), length(nHs)),
              dimnames= list(1:N, 1:6, as.character(rs), as.character(deltas),
                             as.character(nHs)))
  
  for (r in rs) {
    ir = which(rs ==r)
    
    for (delta in deltas) {
      idelta = which(deltas==delta)
    
      for (nH in nHs) {
        inH = which(nH == nHs)
        
        for (i in 1:N) {
          cat("Simulation (", r, ", ", delta, ", ", nH, ") : ", i, "/", N, "\n")
          
          repeat { sim = sim.epi(nH = nH, sample.n=nH/10, mu=0, vcst=0, r=r, delta=delta*r); if(sim$sampled.time > 0) break; }
          res[i, , ir, idelta, inH] = TajimaD(sim)
        }
      }
    }
  }
  
  res
}

# alignment generator for simulated whbh
whbh2dna = function(whbhsim, genome.size=10) {

  # list of unique strains obeserved during the simulation 
  strains = unique(whbhsim$sampled.strain)
  
  # Seedy structures
  librstrains = c()
  libr = list()
  nuc = list()
  newsite = 0
  
  # add root of the tree
  libr[[1]] = c(NA)
  nuc[[1]] = c(NA)
  librstrains = c(1)
  
  # generate fake genome
  genome = sample.int(4, genome.size, replace=T)
  
  # iterate on strains to make up genetic information
  for (strain in strains) {
    
    # get path to root
    path = c(strain)
    node = strain
    while(whbhsim$mutation.tree[node] != 0) {
      node = whbhsim$mutation.tree[node]
      path = c(node,path)
    }
    
    # build genetic information of path
    # from root ot leaf
    for (node in path) {
      
      if (node %in% librstrains) next;
      
      # add node to librstrains
      librstrains = c(librstrains, node)
      
      # get its index
      inode = length(librstrains)
      
      # add a new seggregating site
      newsite = newsite + 1
      
      # if we have a dad and it is not the root get the mutations from the dad
      dad = whbhsim$mutation.tree[node]
      
      if (dad != 1) {
        idad = which(librstrains==dad)
        libr[[inode]] = c(libr[[idad]], newsite)
        nuc[[inode]] = c(nuc[[idad]], sample(setdiff(1:4, genome[newsite]),1))
      } else {
        libr[[inode]] = newsite
        nuc[[inode]] = sample(setdiff(1:4, genome[newsite]),1)
      }
    }
  }
  # now librstrains, libr, nuc contain all the information to build a seedy object
  # reduce libary to the subset we are interested in
  #libr = libr[librstrains]
  librtoDNA(whbhsim$sampled.strain, libr, nuc, genome, librstrains, filename="testdna", format="fasta")
}





# pop. gen. statistics for seedy's Wright-Fisher simulations

W = function(s) {
  
  all.sites = unlist(s$libr)
  
  # remove NA form the wild
  all.sites= all.sites[!is.na(all.sites)]
  
  # frequencies of sites
  sites.freqs = table(all.sites)
  
  # find which sites appear in all genotypes
  common = as.numeric(names(sites.freqs)[sites.freqs == length(s$libr)])
  
  # remove common sites fron list of unique sites
  seg.sites = setdiff(unique(all.sites), common)
  
  # number of segegrating sitea
  S = length(seg.sites)

  # sample size
  N = length(s$obs.strain)
  
  # Nei's nucleotide diversity
  pi = meansnps(s$obs.strain, rep(1, N), libr=s$libr, nuc=s$nuc, key=s$librstrains)
  # implementation by Seedy is bogus: frequencies are used but sum is done on
  # lower triangle
  pi = 2*N/(N-1) * pi
  
  # tajima's D
  a1 = sum(1/(1:(N-1)))
  a2 = sum(1/(1:(N-1))^2)
  b1 = (N+1)/3/(N-1)
  b2 = 2 * (N^2 + N + 3)/9/N/(N-1)
  c1 = b1 - 1/a1
  c2 = b2 - (N+2)/a1/N + a2/a1^2
  e1 = c1 / a1
  e2 = c2 / (a1^2 + a2)
  
  # Watterson estimator
  W = S/a1 
  
  # Tajima's D
  D = ifelse(S==0, 0, (pi - W)/sqrt(e1 * S + e2 * S * (S-1)))
  
  return(c(pi=pi, W=W, D=D))
}

W.seedy.epi = function(epi,sample.time) {
  dat = epi$sampledata
  s = list(obs.strain=dat[dat[,"sampletimes"]==sample.time,"sampleWGS"],
           libr = epi$libr,
           nuc = epi$nuc,
           librstrains = epi$librstrains)
  W(s)
}
# ##### Wright Fisher simulations
# 
# 
# # plots for Mark
# par(mfrow=c(3,2))
# # Watterson
# plot(density(Ws[2,]), main="WF Watterson")
# plot(density(res.all[,2,1,2,2], na.rm=T), main="WHBH Watterson")
# 
# # Nei's nucleotide diversity
# plot(density(Ws[1,]), main="WF: Nei's pi")
# plot(density(res.all[,1,1,2,2], na.rm=T), main="WHBH: Nei's pi")
# 
# # Tajima's D
# plot(density(Ws[3,]), main="WF: Tajima's D")
# plot(density(res.all[,3,1,2,2], na.rm=T), main="WHBH: Tajima's D")
# 
# #plots different parameters
# rs= c(1, 1.5, 2)
# for (ir in 1:3) {
# # Watterson
# plot(density(res.all[,2,ir,2,1], na.rm=T), main="WHBH Watterson peak")
# plot(density(res.all[,5,ir,2,1], na.rm=T), main="WHBH Watterson")
# 
# # Nei's nucleotide diversity
# plot(density(res.all[,1,ir,2,1], na.rm=T), main="WHBH: Nei's pi - peak")
# plot(density(res.all[,4,ir,2,1], na.rm=T), main="WHBH: Nei's pi")
# 
# # Tajima's D
# 
# plot(density(res.all[,3,ir,2,1], na.rm=T), main="WHBH: Tajima's D - peak")
# plot(density(res.all[,6,ir,2,1], na.rm=T), main="WHBH: Tajima's D")
# }
# 
# # a quick look at corrected values of pi
# #
# res.cor = res.all[,1:3,3,2,2]
# stdev.D = (res.cor[,1]-res.cor[,2])/res.cor[,3]
# hist((res.cor[,1] * 999/1000 - res.cor[,2])/stdev.D)
# hist(res.cor[,3])


# debug of TajimaD
#s = sim.epi(1000,sample.n=1000, mu=0, vcst=0, mut=1e-2)
#whbh2dna(s,100)
#tajima.test(read.FASTA("testdna"))
#TajimaD(s)