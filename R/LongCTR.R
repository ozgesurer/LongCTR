LongCTR <- function(X, L, y, Kmax = 60, cv = 10, warm.start = NULL, family = NULL){
  groupobj <- function(id = id, super_id = super_id, pred_id = c(1:p),
                       lb = seq(1, L*p, by = L), ub = seq(L, L*p, by = L),
                       t_lb = 1, t_ub = L, is_included = FALSE){
    g <- list(id = id, super_id = super_id, pred_id = pred_id, lb = lb, ub = ub, t_lb = t_lb, t_ub = t_ub,
              size = length(pred_id), period = t_ub - t_lb + 1, is_included = is_included)
    class(g) <- "group"
    return(g)
  }
  grouplistobs <- function(p, L){
    gs = list(groupobj(id = 1, super_id = 0, pred_id = c(1:p), lb = seq(1, L*p, by = L), ub = seq(L, L*p, by = L),
                       t_lb = 1, t_ub = L, is_included = FALSE))
    class(gs) = "grouplist"
    return(gs)
  }
  groupstrobj <- function(X, y, L, K){

    n = dim(X)[1]
    p = dim(X)[2]/L
    X = X
    y = y
    A = NA
    K = K
    Z = matrix(1, nrow = n, ncol = 1)
    basis = solve(t(Z)%*%Z)

    E = X - Z %*% (basis %*% (t(Z) %*% X))
    u = colSums(E*E)
    w = t(E)%*%y
    r = ifelse(u == 0, 0, sign(w)*(w^2/u))

    size = 1
    alpha = 0

    SSE = sum((y - mean(y))^2)
    grstr <- list(n = n, p = p, L = L, X = X, y = y, E = E, u = u, w = w, r = r, Z = Z, A = A, K = K,
                  alpha = alpha, alphaorg = 0, basis = basis, groups = grouplistobs(p, L), size = size, SSE = SSE)

    class(grstr) = "groupstr"
    return(grstr)
  }
  grouplist.update <- function(groups, size, new_groups, X){

    n <- dim(X)[1]
    total_p <- dim(X)[2]

    super_id <- new_groups$split_group$id

    for (i in 1:length(groups)){
      groups[[i]]$super_id = groups[[i]]$id
    }
    g1 = groupobj(id = super_id, super_id = super_id, pred_id = new_groups$pred_in,
                  lb = new_groups$mid, ub = new_groups$ub,
                  t_lb = new_groups$t_mid, t_ub = new_groups$t_ub,
                  is_included = TRUE)
    groups[[g1$id]] <- g1
    size <- size + 1
    z1 <- rowSums(matrix(X[, rep(g1$lb, each = g1$period) + 0:(g1$period - 1)], nrow = n))
    a1 <- rep(0, total_p)
    a1[rep(g1$lb, each = g1$period) + 0:(g1$period - 1)] <- 1
    derivedpredictors <- list(z1 = z1, a1 = a1)

    if(new_groups$pred_out[1] != "NA"){
      is_included = ifelse(new_groups$split_group$is_included == TRUE, TRUE, FALSE)
      g2 = groupobj(id = size, super_id = super_id, pred_id = new_groups$pred_out,
                    lb = new_groups$mid_out, ub = new_groups$ub_out,
                    t_lb = new_groups$t_mid, t_ub = new_groups$t_ub,
                    is_included = is_included)
      groups[[size]] = g2
      size = size + 1
      if(is_included == TRUE & (g1$period != new_groups$split_group$period)){
        z2 <- rowSums(matrix(X[, rep(g2$lb, each = g2$period) + 0:(g2$period - 1)], nrow = n))
        a2 <- rep(0, total_p)
        a2[rep(g2$lb, each = g2$period) + 0:(g2$period - 1)] <- 1
        derivedpredictors <- list(z1 = z1, z2 = z2, a1 = a1, a2 = a2)
      }
      if(g1$period != new_groups$split_group$period){
        is_included = ifelse(new_groups$split_group$is_included == TRUE, TRUE, FALSE)
        g3 = groupobj(id = size, super_id = super_id, pred_id = c(new_groups$pred_in, new_groups$pred_out) ,
                      lb = c(new_groups$lb, new_groups$lb_out), ub = c(new_groups$mid - 1, new_groups$mid_out - 1),
                      t_lb = new_groups$t_lb, t_ub = new_groups$t_mid - 1,
                      is_included = is_included)
        groups[[size]] = g3
      }
    }else{
      if(g1$period != new_groups$split_group$period){
        is_included = ifelse(new_groups$split_group$is_included == TRUE, TRUE, FALSE)
        g3 = groupobj(id = size, super_id = super_id, pred_id = new_groups$pred_in,
                      lb = new_groups$lb, ub = new_groups$mid - 1,
                      t_lb = new_groups$t_lb, t_ub = new_groups$t_mid - 1,
                      is_included = is_included)
        groups[[size]] = g3
      }
    }
    return(list(groups, derivedpredictors))
  }
  group.split <- function(group, E, w, k, n){
    size = group$size
    period = group$period
    t_ub = group$t_ub
    t_lb = group$t_lb

    RSums = rep(0, size)
    vSums = rep(0, size)
    wSums = rep(0, size)
    ESums = matrix(0, nrow = n, ncol = size)

    Rbest = 0

    for(t in 1:period){
      for(j in 1:size){
        pred_lb <- group$ub[j] - t + 1
        pred_ub <- group$ub[j]

        vSums[j] <- sum(w[c(pred_lb:pred_ub)])
        ez <- rowSums(matrix(E[, (pred_lb:pred_ub)], nrow = n))
        ESums[, j] <- ez
        wSums[j] <- t(ez)%*%ez
        RSums[j] <- ifelse(vSums[j]^2 > 1e-10 & wSums[j] > 1e-10, sign(vSums[j])*vSums[j]^2/wSums[j], 0)
      }
      orderedsums <- order(RSums, decreasing = TRUE)
      if(abs(RSums[orderedsums[1]]) < abs(RSums[orderedsums[size]]))
        orderedsums <- rev(orderedsums)

      if(size == 1){
        cumsumE = ESums
      }else{
        cumsumE = t(apply(ESums[, orderedsums], 1, cumsum))
      }
      numerator = cumsum(vSums[orderedsums])
      denominator = colSums(cumsumE * cumsumE)
      R = ifelse(((numerator^2 > 1e-10) & (denominator > 1e-10)), (numerator^2)/denominator, 0)
      row_star = which.max(R)

      if(R[row_star] > Rbest){
        Rbest = R[row_star]
        col_star = t

        pred_in_gr = group$pred_id[orderedsums][1:row_star]
        ub_in_gr = group$ub[orderedsums][1:row_star]
        lb_in_gr = group$lb[orderedsums][1:row_star]

        if (row_star == size){
          pred_not_in_gr = "NA"
          ub_out = 0
          lb_out = 0
        }else{
          pred_not_in_gr = group$pred_id[orderedsums][(row_star + 1):size]
          ub_out = group$ub[orderedsums][(row_star + 1):size]
          lb_out = group$lb[orderedsums][(row_star + 1):size]
        }
      }
    }

    if(Rbest > 0){
      return(list(split_group = group, red = Rbest, pred_in = pred_in_gr, pred_out = pred_not_in_gr,
                  lb = lb_in_gr, mid = ub_in_gr - col_star + 1, ub = ub_in_gr,
                  lb_out = lb_out, mid_out = ub_out - col_star + 1, ub_out = ub_out,
                  t_lb = t_lb, t_mid = t_ub - col_star + 1, t_ub = t_ub))
    }else{
      return(list(split_group = group, red = Rbest))
    }
  }
  groupstr.update <- function(z, a, groupstr, groups){
    groupstr$groups <- groups
    groupstr$size <- length(groups)
    Z = groupstr$Z
    no_dp = dim(Z)[2]
    basis = groupstr$basis
    y = groupstr$y
    E = groupstr$E
    w = groupstr$w
    u = groupstr$u
    A = groupstr$A

    ez <- z - (Z %*% basis) %*% (t(Z) %*% z)
    ezezinv <- as.numeric(1/(t(ez) %*% ez))
    Z <- cbind(Z, z)
    if(no_dp == 1)
      A <- a
    else
      A <- cbind(A, a)
    basis <- solve(t(Z)%*%Z)
    Zy <- t(Z) %*% y
    alpha <- basis%*%Zy
    Yhat = Z %*% alpha
    SSE = sum((y - Yhat)^2)

    groupstr$SSE = SSE
    v <- t(E) %*% ez
    groupstr$E <- E - ezezinv * (ez %*% t(v))
    groupstr$w <- w - v * as.numeric((t(ez)%*%y) * ezezinv)
    groupstr$u <- u - ezezinv * (v * v)
    groupstr$Z <- Z
    groupstr$A <- A
    groupstr$basis <- basis
    groupstr$alphaorg <- alpha
    groupstr$alpha <- compute.alpha(alpha, groupstr)
    return(groupstr)
  }
  compute.alpha <- function(alpha, groupstr){
    L = groupstr$L
    p = groupstr$p
    totalp = L*p
    A = matrix(groupstr$A, nrow = totalp)
    beta = matrix(c(alpha[1], A %*% alpha[-1]), nrow = totalp + 1)

    groups = groupstr$groups
    alpha = rep(0, length(groups) + 1)
    alpha[0] = beta[0]
    for(j in 1:length(groups)){
      for(i in 1:totalp){
        for(l in 1:length(groups[[j]]$lb)){
          if((groups[[j]]$lb[l] <= i) & (i <= groups[[j]]$ub[l])){
            alpha[j + 1] <- beta[i + 1]
          }
        }
      }
    }
    return(alpha)
  }
  CTRgroup <- function(X, y, groupstr){
    E = groupstr$E
    r = groupstr$r
    u = groupstr$u
    w = groupstr$w
    groups = groupstr$groups
    size = groupstr$size
    Z = groupstr$Z
    basis = groupstr$basis
    n = groupstr$n
    p = groupstr$p
    Rbest = 0

    splits_best = list(red = 0)
    for (k in 1:size){
      group <- groups[[k]]
      new_groups <- group.split(group, E, w, k, n)
      if(new_groups$red > Rbest){
        splits_best <- new_groups
        Rbest <- new_groups$red
      }
    }

    if(splits_best$red > 0){
      gsu <- grouplist.update(groups, size, splits_best, X)
      groups <- gsu[[1]]
      derived <- gsu[[2]]
      z <- derived$z1
      a <- derived$a1
      groupstr <- groupstr.update(z, a, groupstr, groups)
      if(length(derived) > 2){
        z <- derived$z2
        a <- derived$a2
        groupstr <- groupstr.update(z, a, groupstr, groups)
      }
    }else{

    }
    return(groupstr)
  }
  groupstr.add <- function(groupstr, groupstrlist, iteration){
    groupstrlist[[iteration]] <- groupstr
    return(groupstrlist)
  }
  if(cv > 0){
    #Run 10-fold CV
    istest <- TRUE
    SSE_cv <- matrix(0, nrow = cv, ncol = Kmax)
    SSE_cv_test <- matrix(0, nrow = cv, ncol = Kmax)
    foldid <- sample(cut(seq(1,  dim(X)[1]), breaks = cv, labels = FALSE));
    for(ifold in 1:cv){
      test <- which(foldid == ifold, arr.ind = TRUE)
      train <- which(foldid != ifold, arr.ind = TRUE)
      Xtest <- X[test, ]
      ytest <- y[test]
      Xtrain <- X[train, ]
      ytrain <- y[train]
      groupstr <- groupstrobj(Xtrain, ytrain, L, Kmax)
      SST <- groupstr$SSE
      for (iteration in 1:Kmax){
        groupstr <- CTRgroup(Xtrain, ytrain, groupstr)
        SSE_cv[ifold, iteration] <- 1 - groupstr$SSE/SST
        yhat <- predict.groupstr(groupstr, Xtest)
        SSE_cv_test[ifold, iteration] <- sum((ytest - yhat)^2)
      }
    }
    # print(colSums(SSE_cv_test))
    K <- which.min(colSums(SSE_cv_test))[1]
    groupstr <- groupstrobj(X, y, L, K)
    groupstrlist <- list()
    groupstrlist <- groupstr.add(groupstr, groupstrlist, 1)
    for (iteration in 1:K){
      groupstr <- CTRgroup(X, y, groupstr)
      groupstrlist <- groupstr.add(groupstr, groupstrlist, iteration + 1)
    }
  }else{
    #Run CTR without CV
    K <- Kmax
    SSE_cv <- rep(0, K)
    groupstr <- groupstrobj(X, y, L, K)
    groupstrlist <- list()
    groupstrlist <- groupstr.add(groupstr, groupstrlist, 1)
    SST <- groupstr$SSE
    for (iteration in 1:K){
      groupstr <- CTRgroup(X, y, groupstr)
      SSE_cv[iteration] <- groupstr$SSE
      groupstrlist <- groupstr.add(groupstr, groupstrlist, iteration + 1)
    }
  }
  return(list(groupstr = groupstr, groupstrlist = groupstrlist, SSE_cv = SSE_cv))
}

predict.groupstr <- function(groupstr, newdata){
  L = groupstr$L
  p = groupstr$p
  totalp = L*p
  A = matrix(groupstr$A, nrow = totalp)
  beta = matrix(c(groupstr$alphaorg[1], A %*% groupstr$alphaorg[-1]), nrow = totalp + 1)
  yhat = cbind(1, newdata) %*% beta
  return(yhat)
}
print.grpstr <- function(groupstrlist, L, K){
  cat("(Group)(alpha){PredIds}{TimeIds}")
  cat("\n")
  #cat( "(0.1)","{", 1, "-", p, "}", sep = "" )
  #cat( "{", 1, "-", L, "}", sep = "" )
  #cat("\n")
  plist <- list()
  totalsum <- 0
  for(k in 1:K){
    groups <- groupstrlist[[k]]$groups
    totalsum <- totalsum + length(groups)
    plist[k] <- list(rep(0, length(groups)))
  }

  print.level <- function(k, j, group, groupstrlist){
    cat(character(k + 1), collapse = " ")
    cat( "(",k, ".", j, ")", "(", round(groupstrlist[[k]]$alpha[j + 1], 4), ")", "{", sep = "")
    cat(sort(group$pred_id), sep = "," )
    cat( "}", sep = "")
    cat( "{", group$t_lb, "-", group$t_ub, "}", sep = "" )
    cat("\n")

    if(k + 1 > K){
      super_id = group$super_id
    }else{
      super_id = group$id
    }
    return(super_id)
  }
  k <- 1
  groups <- groupstrlist[[k]]$groups
  j <- 1
  super_id <- 0
  id <- 1
  terminate = F
  while(terminate == F){
    if(plist[[k]][j] == 0){
      while(plist[[k]][j] == 0){
        if (groups[[j]]$super_id == super_id){
          super_id = print.level(k, j, groups[[j]], groupstrlist)
          plist[[k]][j] = 1
          k <- k + 1
          if(k > K){
            k <- k - 1
            for (t in 1:length(plist[[k]])){

              if((plist[[k]][t] == 0) & (groups[[t]]$super_id == super_id)){
                super_id = print.level(k, t, groups[[t]], groupstrlist)
                plist[[k]][t] = 1
              }
            }

            gg = FALSE
            while (gg == FALSE & k - 1 > 0){
              k <- k - 1
              groups <- groupstrlist[[k]]$groups
              super_id <- groups[[super_id]]$super_id
              for (t in 1:length(plist[[k]])){
                if((plist[[k]][t] == 0) & (groups[[t]]$super_id == super_id)){
                  gg = TRUE
                }
              }
            }
          }
          groups <- groupstrlist[[k]]$groups
          j <- 1
        }else{
          j <- j + 1
          if(j > length(plist[[k]])){
            k <- 1
            j <- 1
            groups <- groupstrlist[[k]]$groups
            super_id <- groups[[j]]$super_id
            if(sum(unlist(plist)) == totalsum){
              terminate = T
            }else{
              while(sum(plist[[k]]) == length(plist[[k]])){
                k <- k + 1
                j <- which(plist[[k]] == 0)[1]
                groups <- groupstrlist[[k]]$groups
                super_id <- groups[[j]]$super_id
              }
            }
          }
        }
      }
    }else{
      j <- j + 1
      if(j > length(plist[[k]])){
        k <- 1
        j <- 1
        groups <- groupstrlist[[k]]$groups
        super_id <- groups[[j]]$super_id
        if(sum(unlist(plist)) == totalsum){
          terminate = T
        }else{
          while(sum(plist[[k]]) == length(plist[[k]])){
            k <- k + 1
            j <- which(plist[[k]] == 0)[1]
            groups <- groupstrlist[[k]]$groups
            super_id <- groups[[j]]$super_id
          }
        }
      }
    }
  }
}
summary.grpstr <- function(groupstrlist,...){

  cat("Coefficient Tree Regression\n\n")
  groups <- groupstrlist[[length(groupstrlist)]]$groups
  cat("Final group structure:\n\n")
  cat("(Group)", "(alpha)", "(Pred Ids)", "(Time Ids)\n", sep = "")


  for (j in 1:length(groups)){
    group = groups[[j]]
    cat( "(", j, ")", "(",round(groupstrlist[[length(groupstrlist)]]$alpha[j + 1], 4), ")", "{", sep = "")
    cat(sort(group$pred_id), sep = "," )
    cat( "}", sep = "")
    cat( "{", group$t_lb, "-", group$t_ub, "}", sep = "" )
    cat( "(", group$is_included,  ")", sep = "" )
    cat("\n")
  }

}
print.tree <- function(groupstrlist, K){

  size.tree <- 0
  for(k in 1:K){
    groups <- groupstrlist[[k]]$groups
    size.tree <- size.tree + length(groups)
  }

  p_lab <- 1
  c_lab <- 2
  c_lab_list <- list('{1, 2, 3}, [1, 10]')
  fromlist <- list()
  tolist <- list()
  alpha_list <- list()

  for(k in 1:(K-1)){
    groups_p <- groupstrlist[[k]]$groups
    groups_c <- groupstrlist[[k+1]]$groups
    for (g in 1:length(groups_p)){
      p_id <- groups_p[g][[1]]$id
      for (g_c in 1:length(groups_c)){
        c_sid <- groups_c[g_c][[1]]$super_id
        if (c_sid == p_id){
          fromlist = append(fromlist, p_lab)
          tolist = append(tolist, c_lab)

          c_lab_list = append(c_lab_list,
                              paste('{', paste(as.character(groups_c[g_c][[1]]$pred_id), collapse=','), '}',
                                    '[', as.character(groups_c[g_c][[1]]$t_lb), ',',
                                    as.character(groups_c[g_c][[1]]$t_ub), ']', collapse=','))

          alpha_list <- append(alpha_list, groupstrlist[[k+1]]$alpha[g_c + 1])

          c_lab = c_lab + 1
        }
      }
      p_lab = p_lab + 1
    }
  }

 #for (k in 2:K){
  #  alpha_list <- append(alpha_list, groupstrlist[[k]]$alpha[2:length(groupstrlist[[k]]$alpha)])
  #}
  ndf <- create_node_df(
    n = size.tree,
    label = unlist(c_lab_list),
    fixedsize = F,
    shape = "rectangle")
  ndf


  edf <- create_edge_df(
    from = unlist(fromlist),
    to   = unlist(tolist),
    rel = c("leading_to"),
    label = paste(round(unlist(alpha_list),2), sep=" "))
  edf

  the_graph <- create_graph(
    nodes_df = ndf,
    edges_df = edf, directed = TRUE)

  render_graph(graph = the_graph, layout = "tree", output = "graph")
}

