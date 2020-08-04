library(gtools)

inclusion_exclusion <- function(event_probs, num_events) {
  probs = c()
  total_events = length(event_probs)
  for (i in 1:total_events) {
    event_combos = combn(event_probs, i, simplify = FALSE)
    for (event_combo in event_combos) {
      probs = c(probs, ((1 - sum(event_combo))**num_events) * ((-1) ** ( (i - 1) %% 2)) )
    }
  }
  return(probs)
}


log_prob_observed <- function(pvec, n, eps_pos, eps_neg, obs_vec) {
  p_eq0 = (1 - pvec) ** n
  p_gr0 = 1 - p_eq0
  tp = obs_vec * (1 - eps_pos) * (p_gr0)
  fp = obs_vec * (eps_pos) * (p_eq0)
  tn = (1 - obs_vec) * (1 - eps_neg) * (p_eq0)
  fn = (1 - obs_vec) * (eps_neg) * (p_gr0)
  total = sum(tp + fp + tn + fn)
  return(sum(log((tp + fp + tn + fn) / total)))
}


pvec = c(.4, .6)
coi = 3
exp(log_prob_observed(pvec, coi, .05, .15, c(0,0))) +
exp(log_prob_observed(pvec, coi, .05, .15, c(1,0))) +
exp(log_prob_observed(pvec, coi, .05, .15, c(0,1))) +
exp(log_prob_observed(pvec, coi, .05, .15, c(1,1)))



pvec = c(.1, .2, .3, .4)
pvec = rdirichlet(1, rep(1, 10))[1,]
n_trials = 7
x = rmultinom(1, n_trials, pvec)[,1] > 0

res = rmultinom(200000, n_trials, pvec)
res2 = apply(res, 2, function(y) { all((y > 0) == x) })

pred_prob = (1 - sum(inclusion_exclusion(pvec[x]/sum(pvec[x]), n_trials))) * (sum(pvec[x])**n_trials)

mean(res2)
print(pred_prob)


for(i in 1:100) {
  x = rmultinom(1, i, pvec)[,1] > 0
  res = rmultinom(100000, i, pvec)
  res2 = apply(res, 2, function(y) { all((y > 0) == x) })

  pred_prob = (1 - sum(inclusion_exclusion(pvec[x]/sum(pvec[x]), i))) * (sum(pvec[x])**i)
}

