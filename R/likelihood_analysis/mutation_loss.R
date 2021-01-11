intp_geom_range <- seq(0.05, .95, .05)
loss_rate_range <- seq(0.05, .95, .05)
mutation_rate <- 1e-8
max_unobs <- 7

intp <- t(sapply(intp_geom_range, function(x) {
  sapply(1:max_unobs, function(y) { dgeom(y, x) })
}))
intp <- intp / rowSums(intp)





res = list()
for (k in 1:length(loss_rate_range)) {
  res[[k]] = list()
  state_change_matrix <- matrix(c(1 - mutation_rate, mutation_rate,
                                  loss_rate_range[k], 1 - loss_rate_range[k]), ncol = 2, byrow = TRUE)
  for (j in 1:length(intp_geom_range)) {
    tmp <- state_change_matrix
    val <- tmp * intp[j,1]
    for (i in 2:max_unobs) {
      tmp <- tmp %*% state_change_matrix
      val <- val + tmp * intp[j,i]
    }
    res[[k]][[j]] = val
  }
}


lost <- 25 * 22
kept <- 55 * 22

llik <- data.frame(llik = c(), loss_rate = c(), geom_pr = c())
surface <- matrix(nrow=length(loss_rate_range), ncol = length(intp_geom_range))
for (i in 1:length(loss_rate_range)) {
  for (j in 1:length(intp_geom_range)) {
    ll <- c(log(res[[i]][[j]][2,1]) * lost + log(res[[i]][[j]][2,2]) * kept)
    llik <- rbind(llik, data.frame(
      value = ll,
      loss_rate = loss_rate_range[i],
      geom_pr = intp_geom_range[j])
      )
    surface[i, j] = ll
  }
}

# image(x=loss_rate_range, y = intp_geom_range, z = surface)
levelplot(surface, row.values = loss_rate_range, column.values = intp_geom_range, xlim = c(0.05, 0.95), ylim = c(0.05, 0.95), contour = T)

which(intp_geom_range == .5)


qplot(x = loss_rate_range, y = surface[, which(intp_geom_range == .5)])
