#' Function to simulate PPS data 
#' 
#' Simulates PPS data for a prespecified hospital according to a steady state model of
#' incidence ( see vignette "Overview of prevtoinc-Package" for details.)
#' 
#' @param n.sample number of beds to simulate
#' @param steps number of steps to evolve the process
#' @param hospital  type of hospital as a list-object (see vignette for details)
#' @return data frame with following columns \itemize{
#'    \item{A.loi - length of infection up to PPS}
#'    \item{L.loi - total length of infection}
#'    \item{A.los - length of stay up to PPS}
#'    \item{L.los - total length of stay}
#'    \item{patient.type - patient type}
#' }
#'
#'@examples
#'
#' pat.1 <- list(dist.X.los = 
#'                    create_dist_vec(function(x) dpois(x-1, lambda = 12), 70),
#'               I.p = 0.008,
#'               dist.X.loi = 
#'                    create_dist_vec(function(x) dpois(x-1, lambda = 10), 70))
#' 
#' pat.2 <- list(dist.X.los = 
#'                    create_dist_vec(function(x) dpois(x-1, lambda = 10), 70),
#'               I.p = 0.02,
#'               dist.X.loi =
#'                    create_dist_vec(function(x) dpois(x-1, lambda = 7), 70))
#' 
#' patient.list <- list(pat.1, pat.2)
#' 
#' 
#' # define distribution of patients
#' pat.1.prob <- 0.4; pat.2.prob <- 0.6
#' pat.dist.hosp <- c(pat.1.prob, pat.2.prob)
#' hospital.1 <- list(inc.factor = 1,
#'                    pat.dist = pat.dist.hosp,
#'                    patient.list = patient.list)
#' data.pps <- simulate_pps_data(n.sample=1000, steps=200, hospital=hospital.1)                    
#'
#' @export
simulate_pps_data <- function(n.sample, steps, hospital) {
  t(replicate(n.sample, 
              bed.process(steps, hospital))) %>%
    tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(unlist))
}

#' Faster method to generate data for PPS with only length of nosocomial infections as output
#' 
#' The function `simulate_pps_fast` can be used to generate PPS data.
#' This functions simulates a PPS on the basis of a given prevalence `P` using
#' a vector of probabilities `dist.X.loi` for the values 1:length(dist.X.loi) of X.loi.
#' It directly samples the time of infection up to date based on `dist.X.loi`.
#' Optionally, the length of stay is sampled independently ( treating the marginal
#' distributions of length of stay and length of infection as independent
#' by assumption) using `dist.X.los` which is in the same format as `dist.X.loi`. 
#' Because of this non-joint sampling rows should not be interpreted as individual 
#' patients.
#' 
#' @param n.sample number of beds to simulate
#' @param P prevalence of nosocomial infections
#' @param dist.X.loi vector of probabilities for values 1:length(dist.X.loi) of X.loi
#' @param dist.X.los vector of probabilities for values 1:length(dist.X.los) of X.los
#' @param one.factor.loi factor by which to approx. reduce number of one day observations for A.loi
#' @param one.factor.los factor by which to approx. reduce number of one day observations for A.los
#' @return data frame with a row for a each simulated patient and  the following columns 
#' \itemize{ 
#'    \item{A.loi - length of infection up to PPS}
#'    \item{L.loi - total length of infection}  
#'    \item{A.los - length of stay up to PPS} 
#'    \item{L.los - total length of stay} 
#'    \item{patient.type - patient type (fixed to 1 for fast method)}
#' }
#' @examples
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=5000,
#'                                   P=0.05,
#'                                   dist.X.loi = example.dist)
#' head(data.pps.fast)
#' @export
simulate_pps_fast <- function(n.sample,
                              P,
                              dist.X.loi,
                              dist.X.los = NA,
                              one.factor.loi = 1,
                              one.factor.los = 1) {
  # simulate A.loi and L.loi
  max.dist.loi <- length(dist.X.loi)
  n.noso <- rbinom(1,n.sample, P)
  x.loi <- sum(1:max.dist.loi*dist.X.loi)
  p.L.loi <- dist.X.loi*(1:max.dist.loi)/x.loi
  L.loi <- sample(1:max.dist.loi, n.noso, T, p.L.loi)
  A.loi <- sapply(L.loi, function(i) sample.int(i,1))
  
  # deflate 1 day A.loi events by factor
  if ( one.factor.loi < 1 ) {
    nbr.ones <- sum(A.loi==1)
    ones.to.cut <- floor((1-one.factor.loi)*nbr.ones)
    A.loi <- sort(A.loi)
    A.loi <- c(rep(0,ones.to.cut), A.loi[(ones.to.cut+1):length(A.loi)])
  }
  A.loi <- c(rep(0, n.sample-n.noso), A.loi)
  L.loi <- c(rep(0, n.sample-n.noso), L.loi)
  
  # simulate A.los and L.los (optionally)
  if (any(!is.na(dist.X.los))) {
    # simulate A.los and L.los
    max.dist.los <- length(dist.X.los)
    x.los <- sum(1:max.dist.los*dist.X.los)
    p.L.los <- dist.X.los*1:max.dist.los/x.los
    L.los <- sample(1:max.dist.los, n.sample, T, p.L.los)
    A.los <- sapply(L.los, function(i) sample.int(i,1))
    
    # deflate 1 day A.los events by factor
    if ( one.factor.los < 1 ) {
      nbr.ones <- sum(A.los==1)
      ones.to.cut <- floor((1-one.factor.los)*nbr.ones)
      A.los <- sort(A.los)
      L.los <- sort(L.los)
      A.los <- A.los[(ones.to.cut+1):length(A.los)]
      L.los <- L.los[(ones.to.cut+1):length(L.los)]
      L.los.add <- sample(2:max.dist.los, 
                          ones.to.cut,
                          T, 
                          p.L.los[2:max.dist.los]/
                            sum(p.L.los[2:max.dist.los]))
      A.los.add <- sapply(L.los.add-1, function(i) sample.int(i,1))+1
      A.los <- c(A.los, A.los.add)
      L.los <- c(L.los, L.los.add)
    }
  } else {
    A.los <- NA
    L.los <- NA
  }
  
  dplyr::data_frame(
    A.loi,
    L.loi,
    A.los,
    L.los,
    patient.type = 1)
}

#' Calculate theoretical values like 
#' x.los, x.loi and other characteristics of the patient population
#' 
#'
#' @param hospital type of hospital as a list-object (see vignette for details)
#' @param steps number of steps to evolve process
#'
#' @return list with following components \itemize{
#'    \item{x.los - average length of stay x_{los}} 
#'    \item{x.loi - average length of infection x_{loi}}
#'    \item{x.los.noso.only - average length of stay for patients with HAI}
#'    \item{x.los.wo.noso - average length of stay for patients discounting time with HAI}    
#'    \item{I -  theoretical incidence rate per patient day}
#'    \item{I.pp - list of theoretical incidences for patient types}
#'    \item{patient.stats - list with `x.los` and `x.loi` for different patient types}
#'    \item{patient.risk.times - list of patient days at risk for different patient types}
#' }
#' @examples 
#' pat.1 <- list(dist.X.los = create_dist_vec(
#'                                function(x) dpois(x-1, lambda = 12), 70),
#' I.p = 0.008,
#' dist.X.loi = create_dist_vec(function(x) dpois(x-1, lambda = 10), 70))
#' 
#' pat.2 <- list(dist.X.los = 
#'                 create_dist_vec(function(x) dpois(x-1, lambda = 10), 70),
#'               I.p = 0.02,
#'               dist.X.loi = 
#'                 create_dist_vec(function(x) dpois(x-1, lambda = 7), 70))
#' 
#' patient.list <- list(pat.1, pat.2)
#' 
#' 
#' # define distribution of patients
#' pat.1.prob <- 0.4; pat.2.prob <- 0.6
#' pat.dist.hosp <- c(pat.1.prob, pat.2.prob)
#' hospital.1 <- list(inc.factor = 1,
#'                    pat.dist = pat.dist.hosp,
#'                    patient.list = patient.list)
#' data.pps <- simulate_pps_data(n.sample=1000, steps=200, hospital=hospital.1) 
#' data.inc.theo <- simulate_incidence_stats(hospital.1, 365 * 1000)
#' # gives incidence rate I
#' data.inc.theo$I
#' # gives incidence proportion per admission
#' data.inc.theo$I.pp
#' @export
simulate_incidence_stats <- function(hospital, 
                                     steps = 365 * 10000) {
  results <- bed.process.long(steps, hospital)
  x.los <- steps/results$n.patients
  x.loi <- results$noso.state/results$noso.count
  I <- results$noso.count/(steps - results$noso.state)
  x.los.only.noso <- results$noso.to.dist.count/results$n.noso.patients
  patient.risk.times <- results$patient.risk.times
  # length of stay from onset of first noso
  
  patient.list <- hospital$patient.list
  patient.stats <- lapply(1:length(patient.list), function(i) {
    x.loi <- mean(sample(1:length(patient.list[[i]]$dist.X.loi),
                         5000, replace = TRUE,
                         prob = patient.list[[i]]$dist.X.loi))
    x.los <- mean(sample(1:length(patient.list[[i]]$dist.X.los),
                         5000, replace = TRUE,
                         prob = patient.list[[i]]$dist.X.los))
    list(x.loi = x.loi, x.los = x.los)
  })
  
  return(list(x.los = x.los,
              x.loi = x.loi,
              I = I,
              I.pp = results$noso.count/results$n.patients,
              x.los.only.noso = x.los.only.noso,
              x.los.wo.noso = (steps - results$noso.state)/results$n.patients,
              patient.stats = patient.stats,
              patient.risk.times = patient.risk.times)
  )
}

#' Function to calculate theoretical value for x.loi and I
#'
#' @param P prevalence of HAIs 
#' @param dist.X.loi probability mass function of distribution of lengths of infection
#' @param dist.X.los vector of probabilities for values 1:length(dist.X.los) of X.los
#'
#' @return list with following components \itemize{
#'    \item{x.loi - average length of infection}
#'    \item{x.los - average length of stay}
#'    \item{I - theoretical incidence rate per patient day}
#'    \item{I.pp - theoretical incidence proportion per admission}
#' }
#' @examples 
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), 70)
#' data.fast.inc.theo <- simulate_incidence_stats_fast(P=0.05, dist.X.loi = example.dist)
#' data.fast.inc.theo$x.loi
#' data.fast.inc.theo$I
#' @export
simulate_incidence_stats_fast <- function(P, dist.X.loi, dist.X.los = NA) {
  x.loi <-  sum(1:length(dist.X.loi) * dist.X.loi)
  if (any(!is.na(dist.X.los))) {
    x.los <-  sum(1:length(dist.X.los) * dist.X.los)
  } else {
    x.los <- NA
    I.pp <- NA
  }
  list(x.loi = x.loi,
       x.los = x.los,
       I = P/(1-P)/x.loi,
       I.pp = P*x.los/x.loi
  )
}

# function that returns the evolution of a single bed. It
# should return for each time the nights the patient has
# passed so far and the state of noscomial infection 
# Input:
# steps - # of steps to evolve the process 
# hospital - hospital considered 
# Output:
# los - length of stay up to date 
# loi - length of a  possible noso up to date
# noso.total - total length of the nosocomial infection (including time after date)
# patient.type - type of patient occupying bed
bed.process <- function(steps, hospital) {
  patient.list <- hospital$patient.list  # list of possible patients
  pat.dist <- hospital$pat.dist  # distribution of patients for hospital
  inc.factor <- hospital$inc.factor  # hospital-factor for incidence
  A.los <- 0 # counter giving number of days of stay left
  A.loi <- 0 # counter giving number of days of infection left
  
  # starting distribution
  patient.nbr <- sample.int(length(patient.list),
                            1,
                            prob = pat.dist)
  los.time <- sample.int(length(patient.list[[patient.nbr]]$dist.X.los),
                         1,
                         prob = patient.list[[patient.nbr]]$dist.X.los)
  loi.dist.curr <- patient.list[[patient.nbr]]$dist.X.loi
  I.current <- inc.factor * patient.list[[patient.nbr]]$I.p
  L.los <- los.time
  A.los <- 0
  hasHAI <- FALSE
  L.loi <- 0
  A.loi <- 0
  
  
  while (steps > 0) {
    A.los <- A.los + 1
    if (hasHAI) {
      A.loi <- A.loi + 1
      if (A.loi > L.loi) {
        hasHAI <- FALSE
        A.loi <- 0
        L.loi <- 0
      }
    }
    
    if (A.los > L.los) {
      # new patient arrives :
      # base length of stay and patient type determined
      patient.nbr <- sample.int(length(patient.list),
                                1,
                                prob = pat.dist)
      L.los <- sample.int(length(patient.list[[patient.nbr]]$dist.X.los),
                          1,
                          prob = patient.list[[patient.nbr]]$dist.X.los)
      loi.dist.curr <- patient.list[[patient.nbr]]$dist.X.loi
      I.current <- inc.factor * patient.list[[patient.nbr]]$I.p
      A.los <- 1
    }
    
    # check for nosocomial infection
    if ( !hasHAI & rbinom(1, 1, I.current) ) {
      hasHAI <- TRUE
      L.loi <- sample.int(length(patient.list[[patient.nbr]]$dist.X.loi),
                          1,
                          prob = patient.list[[patient.nbr]]$dist.X.loi)
      L.los <- L.los + L.loi
      A.loi <- 1
    }
    steps <- steps - 1
  }
  # return length of stay so far and length of possible
  # nosocomial infection up to stay.
  state <- list(
    A.loi = A.loi,
    A.los = A.los,
    L.loi = L.loi,
    L.los = L.los,
    patient.type = patient.nbr
  )
  return(state)
}

# same as bed process except that we get a history of states
bed.process.long <- function(steps, hospital) {
  patient.list <- hospital$patient.list  # list of possible patients
  patient.risk.times <- rep(0, length(patient.list))
  pat.dist <- hospital$pat.dist  # distribution of patients for hospital
  inc.factor <- hospital$inc.factor  # hospital-factor for incidence
  los.time <- 0; noso.time <- 0; noso.count <- 0
  noso.to.dist.count <- 0
  noso.state <- 0  # counts each step where we have a nosocomial infection
  n.patients <- 0
  n.noso.patients <- 0
  while (steps > 0) {
    if (los.time == 0) {
      # new patient arrives:
      n.patients <- n.patients + 1
      pat.nbr <- sample.int(length(patient.list),
                            1,
                            prob = pat.dist)
      patient.type <- patient.list[[pat.nbr]]
      L.los <- sample.int(length(patient.type$dist.X.los),
                          1,
                          prob = patient.type$dist.X.los)
      los.time <- L.los
      patient.risk.times[pat.nbr] <- patient.risk.times[pat.nbr] +
        los.time
      loi.dist.curr <- patient.type$dist.X.loi
      I.current <- inc.factor * patient.type$I.p
      los <- 0
      noso.total <- 0
      noso.to.dist <- 0
    }
    # check for nosocomial infection
    if (noso.time >= 1) {
      noso.state <- noso.state + 1
      noso.time <- noso.time - 1
    } else if (rbinom(1, 1, I.current)) {
      # introduce variable to count time from first noso to
      # discharge
      if (noso.to.dist == 0) {
        n.noso.patients <- n.noso.patients + 1
        noso.to.dist <- 1
      }
      noso.count <- noso.count + 1
      noso.time <- sample.int(length(patient.type$dist.X.los),
                              1,
                              prob = loi.dist.curr)
      noso.total <- noso.time
    } else {
      noso.total <- 0
      los.time <- los.time - 1
    }
    # print(c(los = los, steps = steps, lo.noso = lo.noso))
    if (noso.to.dist == 1) {
      noso.to.dist.count <- noso.to.dist.count + 1
    }
    los <- los + 1
    steps <- steps - 1
  }
  # return length of stay so far and length of possible
  # nosocomial infection up to stay.
  state <- list(
    A.loi = noso.total - noso.time + 1,
    L.loi = noso.total,
    A.los = los,
    L.los = L.los,
    noso.count = noso.count,
    n.patients = n.patients,
    noso.state = noso.state,
    noso.to.dist.count = noso.to.dist.count,
    n.noso.patients = n.noso.patients,
    patient.risk.times = patient.risk.times
  )
  return(state)
}


# distributions used in paper



#' Probability mass function for a Poisson distribution shifted by one and resulting expected value 8
#' 
#' @param x vector of positive integer values to evaluate
#' 
#' @examples 
#' 
#' plot(pois_dist_fct(1:100))
#' 
#' @export
pois_dist_fct <- function(x) dpois(x-1, 7)


#' Probability mass function for a geometric distribution shifted by one and resulting expected value 8
#' 
#' @param x vector of positive integer values to evaluate
#' 
#' @examples 
#' 
#' plot(geom_dist_fct(1:100)) 
#' 
#' @export
geom_dist_fct <- function(x) dgeom(x-1, 1/8)


#' Create vector of probabilities for a finite positive discrete distribution
#' 
#' Cuts-off the (possibly unbounded) probability distribution at `max.dist`
#' and normalizes the resulting vector of probability to sum up to 1.
#' 
#' @param dist probability mass function to use
#' @param max.dist maximum value at which to cutoff distribution
#'
#' @return vector of probabilites for values 1:max.dist
#' 
#' @examples 
#' 
#'
#' geom_dist_fct(1:70)
#' create_dist_vec(geom_dist_fct, max.dist = 70)
#' 
#' @export
create_dist_vec <- function(dist, max.dist) {
  dist(1:max.dist)/sum(dist(1:max.dist))
}

#' Function to simulate PPS and data and calculate a number of estimators
#' @param n.sample number of beds to simulate
#' @param P  average prevalence of nosocomial infections
#' @param dist.X.loi vector of probabilities for values 1:length(dist.X.loi) of X.loi
#' @param data.theo data frame with theoretical info generated by simulate_incidence_stats_* function
#' @param dist.X.los vector of probabilities for values 1:length(dist.X.los) of X.los
#' @param one.factor.loi factor by which to approx. reduce number of one day observations for A.loi
#' @param one.factor.los factor by which to approx. reduce number of one day observations for A.los
#'
#' @return data frame with following columns \itemize{
#'    \item{n - number of patients sampled} 
#'    \item{n.noso - number of HAIs}
#'    \item{P.hat - estimate of prevalence P}
#'    \item{I.hat - estimate of incidence rate I}
#'    \item{I.pp.hat - estimate of incidence proportion per admission I.pp}
#'    \item{x.loi.hat - estimate of x.loi}
#'    \item{x.los.hat - estimate of x.los}
#'    \item{method - name of the method}
#' }
#' and rows for the estimators gren, rear, pps.median, pps.mean, pps.mixed, 
#' rhame.theo, L.full 
#' (for a description of the estimators see vignette).
#' 
#' 
#' @examples 
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' generate_I_fast(200, P = 0.05, example.dist )
#' 
#' @export
generate_I_fast <- function(n.sample, P,
                            dist.X.loi,
                            data.theo = NULL, 
                            dist.X.los = NA,
                            one.factor.loi = 1,
                            one.factor.los = 1) {
  data <- simulate_pps_fast(n.sample, P, dist.X.loi, 
                            dist.X.los = dist.X.los,
                            one.factor.loi = one.factor.loi,
                            one.factor.los = one.factor.los)
  calculate_I(data, data.theo)
}


#' Function to calculate different estimators for I from PPS data.
#'
#' @param data data frame as generated by `simulate_pps_data` or `simulate_pps_fast`
#' @param data.theo data frame as generated by `simulate_incidence_stats` or `simulate incidence_stats_fast``
#'
#' @return data frame with following columns \itemize{
#'    \item{n - number of patients sampled} 
#'    \item{n.noso - number of HAIs}
#'    \item{P.hat - estimate of prevalence P}
#'    \item{I.hat - estimate of incidence rate I}
#'    \item{I.pp.hat - estimate of incidence proportion per admission I.pp}
#'    \item{x.loi.hat - estimate of x.loi}
#'    \item{x.los.hat - estimate of x.los}
#'    \item{method - name of the method}
#' }
#' and rows for the estimators gren, rear, pps.median, pps.mean, pps.mixed, 
#' rhame.theo, L.full 
#' (for a description of the estimators see vignette).
#' 
#' @examples 
#' 
#' # create example data for PPS
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=200,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' head(data.pps.fast)
#' 
#' # estimate incidence based on Rhame-Sudderth formula
#' calculate_I(data = data.pps.fast)
#' 
#' @export
#'
calculate_I <- function(data, data.theo = NULL) {
  
  if (is.null(data.theo)) {
    data.theo <- list(x.loi = NA, x.los = NA, I = NA, I.pp = NA)
  }
  
  # estimates of x.loi
  A.loi <- na.omit(data$A.loi[data$A.loi > 0])
  x.loi.hat.median <- median(A.loi)
  x.loi.hat.mean <- mean(A.loi)
  if ("L.loi" %in% colnames(data) && any(!is.na(data$L.loi))) {
    L.loi <- na.omit(data$L.loi[data$L.loi > 0])
    x.loi.hat.L.full <- length_unbiased_mean(epmf(L.loi))
  } else {
    x.loi.hat.L.full <- NA
  }
  
  # estimates of x.los
  if ("A.los" %in% colnames(data) && any(!is.na(data$A.los))) {
    A.los <- na.omit(data$A.los[data$A.los > 0])
    x.los.hat.median <- median(A.los)
    x.los.hat.mean <- mean(A.los)
  } else {
    x.los.hat.median <- NA
    x.los.hat.mean <- NA
  }
  if ("L.los" %in% colnames(data) && any(!is.na(data$L.los))) {
    L.los <- na.omit(data$L.los[data$L.los > 0])
    x.los.hat.L.full <- length_unbiased_mean(epmf(L.los))
  } else {
    x.los.hat.L.full <- NA
  }
  
  
  
  I.new.gren  <- calculate_I_smooth(data,
                                    method = "gren")
  I.new.rear  <- calculate_I_smooth(data,
                                    method = "rear")
  I.pps.median <- calculate_I_rhame(data,
                                    x.loi.hat.median,
                                    x.los.hat.median,
                                    method = "pps.median")
  I.pps.mean <- calculate_I_rhame(data,
                                  x.loi.hat.mean,
                                  x.los.hat.mean,
                                  method = "pps.mean")

  
  I.full      <- calculate_I_rhame(data,
                                   x.loi.hat.L.full,
                                   x.los.hat.L.full,
                                   method = "L.full")
  I.rhame     <- calculate_I_rhame(data,
                                   data.theo$x.loi,
                                   data.theo$x.los,
                                   method = "rhame.theo")
  
  I.new.mixed <- calculate_I_mixed(I.pps.mean,
                                   I.new.gren,
                                   a = 0.01,
                                   b = 500)

  # construct mixed estimator
  
  
  dplyr::bind_rows(
    I.new.gren,
    I.new.rear,
    I.new.mixed,
    I.pps.median,
    I.pps.mean,
    I.full,
    I.rhame
  )
}