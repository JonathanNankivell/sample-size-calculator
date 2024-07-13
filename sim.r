# NOTATION
# n: number of patients in trial.
# prior.i: prior number of successes on treatment A.
# prior.j: prior number of failures on treatment A.
# prior.k: prior number of successes on treatment B.
# prior.l: prior number of failures on treatment B.
# t: number of patients that have been treated.
# i: observed number of successes on treatment A.
# j: observed number of failures on treatment A.
# k: observed number of successes on treatment B.
# l: observed number of failures on treatment B.


thomper <- function( i, j, k, l, chooser_params ) {

    sample.c <- rbeta(1, i + prior.i, j + prior.j)
    sample.n <- rbeta(1, k + prior.k, l + prior.l)

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}

balanced <- function( i, j, k, l, chooser_params ) {

    n <- i + j + k + l

    if (n%%2) {
        action <- sample(c(1,2), size=1)
    }
    else {
        if (i + j < k + l) {
            action <- 1
        } else {
            action <- 2
        }
    }

    return(action)

}

greedy_with_inference <- function( i, j, k, l, chooser_params ) {

    expected.c <- (prior.i + i)/(prior.i + i + prior.j + j)
    expected.n <- (prior.k + k)/(prior.k + k + prior.l + l)

    if (expected.c > expected.n) {
        return(1)
    } else {
        return(2)
    }
}

greedy <- function( i, j, k, l, chooser_params ) {

    expected.c <- (prior.i)/(prior.i + prior.j)
    expected.n <- (prior.k)/(prior.k + prior.l)

    if (expected.c > expected.n) {
        return(1)
    } else {
        return(2)
    }
}

fallibilist_thomper <- function( i, j, k, l, chooser_params ) {

    para <- 0.1

    if (sample(c(TRUE,FALSE), size=1, prob=c(para,1-para))) {
        if (sample(c(TRUE,FALSE), size=1)) {
            return(1)
        } else {
            return(2)
        }
    }

    sample.c <- rbeta(1, i + prior.i, j + prior.j)
    sample.n <- rbeta(1, k + prior.k, l + prior.l)

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}

foo <- function(i,j,k,l) {
    # https://gwern.net/doc/statistics/decision/1933-thompson.pdf
    numer <- 0

    for (a in 0:(i)) {
        numer <- numer + (choose(i + k - a, k) * choose(j + l + 1 + a, l))
    }
    denom <- choose(i + j + k + l + 2, k + l + 1)
    sum <- numer / denom

    return(sum)
}

tw_foo <- function( i, j, k, l, chooser_params ) {

    c <- chooser_params

    prob <- tryCatch(
        {
            min(foo(i + prior.i, j + prior.j, k + prior.k, l + prior.l),1)
        },
        error=function(cond) {
            expected.c <- (i + prior.i)/(i + prior.i + j + prior.j)
            expected.n <- (k + prior.k)/(k + prior.k + l + prior.l)

            if (expected.c > expected.n) {
                1
            } else {
                0
            }
        }
    )

    if (is.na(prob)) {
        sample.c <- rbeta(1, i + prior.i, j + prior.j)
        sample.n <- rbeta(1, k + prior.k, l + prior.l)

        if (sample.c > sample.n) {
            prob <- 1
        } else {
            prob <- 0
        }
    }

    quant <- (prob**c)/(prob**c + (1-prob)**c)

    #print(paste("c:",c,"prob:",prob,"quant:",quant,"i:", i, "j:",j ,"k:", k,"l:", l))
    if (is.na(quant)) {
        print(paste(i,j,k,l))
    }

    action <- sample(c(1,2), 1, prob = c(quant,1-quant))

    return(action)
}


# https://github.com/kwathen/ClinicalTrialSimulation/blob/master/Example%202/AnalysisMethods.r
# calculate the probability that one beta dist. is greater than another
IneqCalcBeta <- function(dA1,dB1,dA2,dB2) {
    res <- integrate(fBetaIneqCalc,0,1, dA1 = dA1, dB1 = dB1, dA2 = dA2, dB2 = dB2)
    res$value
}

#Helper functions 
fBetaIneqCalc <- function(x, dA1, dB1, dA2, dB2){x**(dA1-1) * (1-x)**(dB1-1) * pbeta(x,dA2,dB2) / beta(dA1,dB1)}

thomper_interperlated <- function( i, j, k, l, chooser_params ) {
    r <- chooser_params

    run_thomper <- sample(c(1,0), 1, prob = c(r,1-r))

    if (run_thomper) {
        sample.c <- rbeta(1, i + prior.i, j + prior.j)
        sample.n <- rbeta(1, k + prior.k, l + prior.l)

        if (sample.c > sample.n) {
            action <- 1
        } else {
            action <- 2
        }
    } else {
        action <- sample(c(1,2), size=1)
    }

    return(action)
}

thomper_trans <- function( i, j, k, l, chooser_params ) {

    c <- chooser_params[1]
    offset <- chooser_params[2]
    minimum <- chooser_params[3]

    sample.c <- rbeta(1, max(minimum, (i + prior.i)*c + offset), max(minimum, (j + prior.j)*c + offset))
    sample.n <- rbeta(1, max(minimum, (k + prior.k)*c + offset), max(minimum, (l + prior.l)*c + offset))

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}

thomper_falling_min <- function( i, j, k, l, chooser_params ) {

    minimum <- chooser_params[1]
    gradient <- chooser_params[2]
    patient <- i+j+k+l

    sample.c <- rbeta(1, max(minimum - patient*gradient, (i + prior.i)), max(minimum - patient*gradient, (j + prior.j)))
    sample.n <- rbeta(1, max(minimum - patient*gradient, (k + prior.k)), max(minimum - patient*gradient, (l + prior.l)))

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}

ebm <- function( i, j, k, l, chooser_params ) {

    trial_size <- chooser_params
    patient <- i+j+k+l

    if (patient < trial_size) {
        if (patient%%2) {
            return(sample(c(1,2), size=1))
        }
        else {
            if (i + j < k + l) {
                return(1)
            } else {
                return(2)
            }
        }
    } else {
        expected.c <- (prior.i + i)/(prior.i + i + prior.j + j)
        expected.n <- (prior.k + k)/(prior.k + k + prior.l + l)

        if (expected.c > expected.n) {
            return(1)
        } else {
            return(2)
        }
    }
    
}


thomper_log <- function( i, j, k, l, chooser_params ) {

    base <- chooser_params

    sample.c <- rbeta(1, log(i + prior.i, base), log(j + prior.j, base))
    sample.n <- rbeta(1, log(k + prior.k, base), log(l + prior.l, base))

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}

thomper_root <- function( i, j, k, l, chooser_params ) {

    c <- chooser_params

    sample.c <- rbeta(1, (i + prior.i)**c, (j + prior.j)**c)
    sample.n <- rbeta(1, (k + prior.k)**c, (l + prior.l)**c)

    if (sample.c > sample.n) {
        return(1)
    } else {
       return(2)
    }
}


tw <- function( i, j, k, l, chooser_params ) {

    c <- chooser_params

    prob <- tryCatch(
        {
            min(IneqCalcBeta(i + prior.i, j + prior.j, k + prior.k, l + prior.l),1)
        },
        error=function(cond) {
            expected.c <- (i + prior.i)/(i + prior.i + j + prior.j)
            expected.n <- (k + prior.k)/(k + prior.k + l + prior.l)

            if (expected.c > expected.n) {
                1
            } else {
                0
            }
        }
    )

    if (is.na(prob)) {
        sample.c <- rbeta(1, i + prior.i, j + prior.j)
        sample.n <- rbeta(1, k + prior.k, l + prior.l)

        if (sample.c > sample.n) {
            prob <- 1
        } else {
            prob <- 0
        }
    }

    quant <- (prob**c)/(prob**c + (1-prob)**c)

    #print(paste("c:",c,"prob:",prob,"quant:",quant,"i:", i, "j:",j ,"k:", k,"l:", l))
    if (is.na(quant)) {
        print(paste(i,j,k,l))
    }

    action <- sample(c(1,2), 1, prob = c(quant,1-quant))

    return(action)
}


tw_falling <- function( i, j, k, l, chooser_params ) {

    c <- (i+j+k+l+1)/(2*numb_of_patients)

    prob <- tryCatch(
        {
            min(IneqCalcBeta(i + prior.i, j + prior.j, k + prior.k, l + prior.l),1)
        },
        error=function(cond) {
            expected.c <- (i + prior.i)/(i + prior.i + j + prior.j)
            expected.n <- (k + prior.k)/(k + prior.k + l + prior.l)

            if (expected.c > expected.n) {
                1
            } else {
                0
            }
        }
    )

    # prob <- min(IneqCalcBeta(i + prior.i, j + prior.j, k + prior.k, l + prior.l),1)
    quant <- (prob**c)/(prob**c + (1-prob)**c)

    #print(paste("c:",c,"prob:",prob,"quant:",quant,"i:", i, "j:",j ,"k:", k,"l:", l))

    action <- sample(c(1,2), 1, prob = c(quant,1-quant))


    return(action)
}


Neyman_with_updating <- function( i, j, k, l, chooser_params ) {
    expected.c <- (prior.i + i)/(prior.i + i + prior.j + j)
    expected.n <- (prior.k + k)/(prior.k + k + prior.l + l)

    p_neyman <- sqrt(expected.c*(1-expected.c))/(sqrt(expected.c*(1-expected.c))+sqrt(expected.n*(1-expected.n)))

    action <- sample(c(1,2),1,prob=c(p_neyman,1-p_neyman))

    return(action)
}


DTL <- function(i,j,k,l,chooser_params) {
    numb_of_immigrants <- chooser_params
    
    done = FALSE
    while(1-done) {
        urn <- c( 5, 2 + numb_of_immigrants - j, 2 + numb_of_immigrants - l )
        #print(urn)
        ball <- sample(c(0,1,2), 1, prob = c(urn/sum(urn)))
        if (ball==0) {
            numb_of_immigrants <- numb_of_immigrants + 1
        } else {
            done <- TRUE
        }
    }
    return(c(ball,numb_of_immigrants))
}


rosey <- function(i,j,k,l,chooser_params) {
    # Rosenberger et al. (2001a)

    expected.c <- (prior.i + i)/(prior.i + i + prior.j + j)
    expected.n <- (prior.k + k)/(prior.k + k + prior.l + l)

    rosey <- sqrt(expected.c)/(sqrt(expected.c)+sqrt(expected.n))

    action <- sample(c(1,2),1,prob=c(rosey,1-rosey))

    return(action)
}

g <- function(alpha, x, y) {
    if (x==1) {return(0)}
    if (x==0) {return(1)}
    else {
        numer <- y * (y/x)**alpha
        denom <- y * (y/x)**alpha + (1 - y) * ((1 - y)/(1 - x))**alpha
        return(numer/denom)
    }
}

rho <- function(p1, p2) {
    return(((1 - p2)/(1 - p2 + 1 - p1)))
}

DBCD <- function(i,j,k,l,chooser_params) {

    alpha <- chooser_params

    n <- i + j + k + l
    if (n<4) {
        if (n%%2 == 1) {
            return(1)
        } else {
            return(2)
        }
    } else {

        p_hat_m_1 <- (i+0.5)/(i+j+1)
        p_hat_m_2 <- (k+0.5)/(k+l+1)

        q_hat_m_1 <- 1 - p_hat_m_1
        q_hat_m_2 <- 1 - p_hat_m_2

        p_hat_m <- q_hat_m_2 / ( q_hat_m_2 + q_hat_m_1 )

        prob <- g( alpha , (i+j)/n , p_hat_m )

        #print(prob)

        action <- sample(c(1,2),1,prob=c(prob,1-prob))
        return(action)
    }
}


gittin_deterministic <- function( i, j, k, l, chooser_params ) {
    discount_rate <- 0.7
    gittins_index.c <- gittins::bmab_gi_ab(i+prior.i,j+prior.j,gamma=discount_rate,N=100,tol=0.0005,lb=NA,ub=NA)
    gittins_index.n <- gittins::bmab_gi_ab(k+prior.k,l+prior.l,gamma=discount_rate,N=100,tol=0.0005,lb=NA,ub=NA)
    if (gittins_index.c == gittins_index.n) {
        action <- sample(c(1,2),1)
        return(action)
    } else if (gittins_index.c > gittins_index.n) {
        return(1)
    } else {
        return(2)
    }
}


FLGI <- function(i,j,k,l,chooser_params) {

    print(c(i,j,k,l))
    block_size <- chooser_params

    discount_rate <- 0.99
    if (block_size == 2) {
        combinations <- list(c(0,0),c(0,1),c(1,0),c(1,1))
    } else if (block_size == 4) {
       combinations <- list(c(0,0,0,0),c(0,0,0,1),c(0,0,1,0),c(0,0,1,1),
                            c(0,1,0,0),c(0,1,0,1),c(0,1,1,0),c(0,1,1,1),
                            c(1,0,0,0),c(1,0,0,1),c(1,0,1,0),c(1,0,1,1),
                            c(1,1,0,0),c(1,1,0,1),c(1,1,1,0),c(1,1,1,1))
    } else if (block_size == 6) {
       combinations <- list(c(0,0,0,0,0,0),c(0,0,0,0,0,1),c(0,0,0,0,1,0),c(0,0,0,0,1,1),
                            c(0,0,0,1,0,0),c(0,0,0,1,0,1),c(0,0,0,1,1,0),c(0,0,0,1,1,1),
                            c(0,0,1,0,0,0),c(0,0,1,0,0,1),c(0,0,1,0,1,0),c(0,0,1,0,1,1),
                            c(0,0,1,1,0,0),c(0,0,1,1,0,1),c(0,0,1,1,1,0),c(0,0,1,1,1,1),
                            c(0,1,0,0,0,0),c(0,1,0,0,0,1),c(0,1,0,0,1,0),c(0,1,0,0,1,1),
                            c(0,1,0,1,0,0),c(0,1,0,1,0,1),c(0,1,0,1,1,0),c(0,1,0,1,1,1),
                            c(0,1,1,0,0,0),c(0,1,1,0,0,1),c(0,1,1,0,1,0),c(0,1,1,0,1,1),
                            c(0,1,1,1,0,0),c(0,1,1,1,0,1),c(0,1,1,1,1,0),c(0,1,1,1,1,1),
                            c(1,0,0,0,0,0),c(1,0,0,0,0,1),c(1,0,0,0,1,0),c(1,0,0,0,1,1),
                            c(1,0,0,1,0,0),c(1,0,0,1,0,1),c(1,0,0,1,1,0),c(1,0,0,1,1,1),
                            c(1,0,1,0,0,0),c(1,0,1,0,0,1),c(1,0,1,0,1,0),c(1,0,1,0,1,1),
                            c(1,0,1,1,0,0),c(1,0,1,1,0,1),c(1,0,1,1,1,0),c(1,0,1,1,1,1),
                            c(1,1,0,0,0,0),c(1,1,0,0,0,1),c(1,1,0,0,1,0),c(1,1,0,0,1,1),
                            c(1,1,0,1,0,0),c(1,1,0,1,0,1),c(1,1,0,1,1,0),c(1,1,0,1,1,1),
                            c(1,1,1,0,0,0),c(1,1,1,0,0,1),c(1,1,1,0,1,0),c(1,1,1,0,1,1),
                            c(1,1,1,1,0,0),c(1,1,1,1,0,1),c(1,1,1,1,1,0),c(1,1,1,1,1,1))
    } else {Abort()}

    actions <- c()
    prob_of_senario <- c()

    #n <- 2
    #list <- rep(list(0:1), n)
    #grid <- expand.grid(list)
    #for(index in 1:n**2) {
    #    senario <- grid[i,]
    for (senario in combinations) {

        # senario would be something like c(0,1) which would
        # symbolise a failure and then a success

        temp.i <- i
        temp.j <- j
        temp.k <- k
        temp.l <- l

        prob_of_results <- c()

        for (result in senario) {
            #result <- as.double(result)
            #print(result)
            # if result is 1, then we're assuming a success
            # if it's 0, then it's a failure

            # to work out the probabililty of this result we must
            # first work out which action has just been played. This
            # is done using the Gittens library.
            # https://github.com/jedwards24/gittins
            gittins_index.c <- gittins::bmab_gi_ab(temp.i+prior.i,temp.j+prior.j,gamma=discount_rate,N=100,tol=0.0005,lb=NA,ub=NA)
            gittins_index.n <- gittins::bmab_gi_ab(temp.k+prior.k,temp.l+prior.l,gamma=discount_rate,N=100,tol=0.0005,lb=NA,ub=NA)
            # print(paste(temp.i+prior.i,temp.j+prior.j,temp.k+prior.k,temp.l+prior.l,gittins_index.c,gittins_index.n,gittins_index.c < gittins_index.n))
            if (gittins_index.c == gittins_index.n) {
                action <- sample(c(1,2),1)
            } else if (gittins_index.c > gittins_index.n) {
                action <- 1
            } else {
                action <- 2
            }
            actions <- append(actions, action)

            # once we know the action that's been chosen, we work
            # out the probability that we would see the result we're
            # assuming.
            if (action==1) {
                prob <- (temp.i + prior.i)/(temp.i + prior.i + temp.j + prior.j)
            } else {
                prob <- (temp.k + prior.k)/(temp.k + prior.k + temp.l + prior.l)
            }
            # switch if we want P(failure)
            # print(result)
            if (result==0) {
                prob <- 1 - prob
            }
            # update i,j,k,l
            if (action==1 & result==1) {
                temp.i <- temp.i+1
            } else if (action==1 & result==0) {
                temp.j <- temp.j+1
            } else if (action==2 & result==1) {
                temp.k <- temp.k+1
            } else {
                temp.l <- temp.l+1
            }

            prob_of_results <- append(prob_of_results, prob)
        }

        prob_of_senario <- append(prob_of_senario, prod(prob_of_results))
    }

    prob_of_senario <- rep(prob_of_senario, each=block_size)

    prob_of_action_1 <- prob_of_senario * (actions==1)
    prob_of_action_1 <- sum(prob_of_action_1)/block_size

    prob_of_action_1 <- max(0,prob_of_action_1)
    prob_of_action_1 <- min(1,prob_of_action_1)

    # print(length(prob_of_senario))
    # print(length(actions))

    action <- sample(c(1,2),1,prob=c(prob_of_action_1,1-prob_of_action_1))

    return(action)
}


forecasters_runner <- function(i, j, k, l, n, chooser, chooser_params) {

    # avoid side-effects
    i.temp <- i
    j.temp <- j
    k.temp <- k
    l.temp <- l

    # set the true probability of success
    p.c <- rbeta( 1, prior.i + i.temp, prior.j + j.temp )
    p.n <- rbeta( 1, prior.k + k.temp, prior.l + l.temp )

    #print(paste(p.c,p.n))

    actions <- c()

    for (iter in (i.temp+j.temp+k.temp+l.temp):n) {

        # use the chooser function to pick an action
        # in this case, it should return 1 or 2
        action <- chooser( i.temp, j.temp, k.temp, l.temp, chooser_params )
        if (length(action)!=1) {
            chooser_params <- action[2]
            action <- action[1]
        }
        #print(action)

        if (action == 1) {

            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.c)

            if (sample == 1) {
                # success!
                i.temp <- i.temp + 1
            } else if (sample == 0) {
                # failure
                j.temp <- j.temp + 1
            }
        } else {
            
            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.n)

            if (sample == 1) {
                # success!
                k.temp <- k.temp + 1
            } else if (sample == 0) {
                # failure
                l.temp <- l.temp + 1
            }
        }
        actions <- append(actions, action)
    }

    is_right_action["expected.c"] <- (prior.i + i.temp)/(prior.i + i.temp + prior.j + j.temp)
    is_right_action["expected.n"] <- (prior.k + k.temp)/(prior.k + k.temp + prior.l + l.temp)

    #print(paste("i:",i.temp,"j:",j.temp,"k:",k.temp,"l:",l.temp))

    #print(paste(is_right_action["expected.c"], is_right_action["expected.n"]))
    #print(is_right_action["expected.c"] < is_right_action["expected.n"])

    is_appears_control_best <- (is_right_action["expected.c"] > is_right_action["expected.n"])

    return(is_appears_control_best)
}

prop_control_appears_best <- function(i, j, k, l, numb_of_sims, n, chooser, chooser_params) {
    # one_run, chooser

    # prop_appears_control_best <- 0

    # for (iter in 1:numb_of_sims) {
    #     is_appears_control_best <- forecasters_runner( i, j, k, l, n, chooser, chooser_params )
    #     prop_appears_control_best <- prop_appears_control_best + (is_appears_control_best - prop_appears_control_best)/(iter+1) #https://en.wikipedia.org/wiki/Moving_average#Cumulative_average
    # }

    list_is_appears_control_best <- list()
    for (iter in 1:1000) {
        is_appears_control_best <- as.double(forecasters_runner( i, j, k, l, n, chooser, chooser_params ))
        #print(is_appears_control_best)
        list_is_appears_control_best <- append(list_is_appears_control_best, is_appears_control_best)
    }
    
    #print(list_is_appears_control_best)
    prop_appears_control_best <- mean(as.double(list_is_appears_control_best))

    return(prop_appears_control_best)
}

## this function doesn't represent what I hoped it would.
## afaict atm, since the probability that the forecasters give affects the trial that is ran,
## this requires an exhausive dynamic programming solutions. I've not found an alternative.
## Needs more thought I think.
tw_with_forecasters <- function( i, j, k, l, chooser_params ) {

    c <- chooser_params[1]
    n <- chooser_params[2]

    prob <- prop_control_appears_best(i, j, k, l, numb_of_sims, n, tw, c)

    quant <- (prob**c)/(prob**c + (1-prob)**c)

    #print(paste("c:",c,"prob:",prob,"quant:",quant,"i:", i, "j:",j ,"k:", k,"l:", l))
    if (is.na(quant)) {
        print(paste(i,j,k,l))
    }

    action <- sample(c(1,2), 1, prob = c(quant,1-quant))

    return(c(action,c(c,n)))
}





one_run <- function(n, chooser, chooser_params) {

    i <- 0
    j <- 0
    k <- 0
    l <- 0

    # set the true probability of success
    p.c <- rbeta( 1, prior.i, prior.j )
    p.n <- rbeta( 1, prior.k, prior.l )

    #print(paste0("p.c: ", p.c, "p.n: ", p.n))

    actions <- c()

    for (m in 1:n) {

        # use the chooser function to pick an action
        # in this case, it should return 1 or 2
        action <- chooser( i, j, k, l, chooser_params )
        if (length(action)!=1) {
            chooser_params <- action[2]
            action <- action[1]
        }

        if (action == 1) {

            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.c)

            if (sample == 1) {
                # success!
                i <- i + 1
            } else if (sample == 0) {
                # failure
                j <- j + 1
            }
        } else {
            
            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.n)

            if (sample == 1) {
                # success!
                k <- k + 1
            } else if (sample == 0) {
                # failure
                l <- l + 1
            }
        }
        actions <- append(actions, action)
    }
    is_right_action <- ((p.c < p.n) + 1) == actions

    xtab <- as.table(rbind(c(i, j), c(k, l)))
    dimnames(xtab) <- list(
        group = c("grp1", "grp2"),
        success = c("yes", "no")
    )
    #print(xtab)
    null_rejected <- tryCatch(
        {
            test <- prop.test(xtab)
            test["p.value"] < 0.05
        },
        error=function(cond) {
            FALSE
        }
    )
    if (is.na(null_rejected)) {
        null_rejected <- FALSE
    }
    #print(null_rejected)

    is_right_action <- append(is_right_action, null_rejected)

    is_right_action["expected.c"] <- (prior.i + i)/(prior.i + i + prior.j + j)
    is_right_action["expected.n"] <- (prior.k + k)/(prior.k + k + prior.l + l)

    #print(paste("i:",i,"j:",j,"k:",k,"l:",l))

    if (p.c > p.n) {
        is_right_action["s_10"] <- i + j + 0.1*(i + j + k + l)< k + l
        is_right_action["prop_best_treatment"] <- (i+j)/(i+j+k+l)
    } else {
        is_right_action["s_10"] <- i + j > k + l + 0.1*(i + j + k + l)
        is_right_action["prop_best_treatment"] <- (k+l)/(i+j+k+l)
    }

    return(is_right_action)
}





one_big_run <- function(n, chooser, chooser_params) {

    i <- 0
    j <- 0
    k <- 0
    l <- 0

    # set the true probability of success
    p.c <- rbeta( 1, prior.i, prior.j )
    p.n <- rbeta( 1, prior.k, prior.l )

    #print(paste0("p.c: ", p.c, "p.n: ", p.n))

    bins <- 1000

    if (n < bins) {
        bins <- n
    }

    actions <- rep(0, bins)
    counts <- rep(0, bins)
    bin_size <- n / bins

    index <- floor(seq(1,bins+0.9999999,(1/bin_size)))

    #print(cbind(actions,counts))

    
    for (m in 1:n) {

        # use the chooser function to pick an action
        # in this case, it should return 1 or 2
        action <- chooser( i, j, k, l, chooser_params )
        if (length(action)!=1) {
            chooser_params <- action[2]
            action <- action[1]
        }

        if (action == 1) {

            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.c)

            if (sample == 1) {
                # success!
                i <- i + 1
            } else if (sample == 0) {
                # failure
                j <- j + 1
            }
        } else {
            
            # decide whether or not the treatment is a success
            sample <- rbinom(1, 1, p.n)

            if (sample == 1) {
                # success!
                k <- k + 1
            } else if (sample == 0) {
                # failure
                l <- l + 1
            }
        }
        is_right_action <- ((p.c < p.n) + 1) == action
        actions[index[m]] <- actions[index[m]] + is_right_action
        counts[index[m]] <- counts[index[m]] + 1
    }

    actions <- actions/counts

    xtab <- as.table(rbind(c(i, j), c(k, l)))
    dimnames(xtab) <- list(
        group = c("grp1", "grp2"),
        success = c("yes", "no")
    )
    #print(xtab)
    null_rejected <- tryCatch(
        {
            test <- prop.test(xtab)
            test["p.value"] < 0.05
        },
        error=function(cond) {
            FALSE
        }
    )
    if (is.na(null_rejected)) {
        null_rejected <- FALSE
    }
    # power is based on _correctly_ rejecting the null
    # we must ignore rejections where the sign is wrong
    p.c_hat <- i/(i+j)
    p.n_hat <- k/(k+l)
    incorrect_direction <- (p.c_hat > p.n_hat) != (p.c > p.n)
    if (null_rejected & incorrect_direction) {
        null_rejected <- FALSE
    }
    #print(null_rejected)

    actions <- append(actions, null_rejected)

    actions["expected.c"] <- (prior.i + i)/(prior.i + i + prior.j + j)
    actions["expected.n"] <- (prior.k + k)/(prior.k + k + prior.l + l)

    #print(paste("i:",i,"j:",j,"k:",k,"l:",l))

    if (p.c > p.n) {
        actions["s_10"] <- i + j + 0.1*(i + j + k + l)< k + l
        actions["prop_best_treatment"] <- (i+j)/(i+j+k+l)
    } else {
        actions["s_10"] <- i + j > k + l + 0.1*(i + j + k + l)
        actions["prop_best_treatment"] <- (k+l)/(i+j+k+l)
    }

    return(actions)
}





simulator <- function(numb_of_sims, n, runner, chooser, chooser_params) {
    # one_run, chooser

    runs <- array(0, dim=n+5)
    props <- c()

    for (iter in 1:(numb_of_sims)) {
        
        is_right_action <- runner(n, chooser, chooser_params)
        props <- append(props, is_right_action[length(is_right_action)])
        is_right_action <- is_right_action[1:(length(is_right_action))]
        runs <- runs + (is_right_action - runs)/(iter+1)
    }

    return(c(runs, props))
}


simulator_new <- function(numb_of_sims, n, runner, chooser, chooser_params) {
    # one_run, chooser

    ns <- rep(n,numb_of_sims)
    raw <- sapply(ns, FUN=runner, chooser=chooser,chooser_params=chooser_params)
    if (as.character(quote(runner)) == "one_big_run") {
        n <- min(n,1000)
    }
    props <- raw[(n+5),]
    runs <- rowMeans(raw)

    return(c(runs, props))
}



## run sims

prior.i <- 1
prior.j <- 1
prior.k <- 1
prior.l <- 1



numb_of_sims <- 100
numb_of_patients <- 100000

x <- 1:(numb_of_patients)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, greedy, NA)
greedy.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
greedy.results <- sim_data[1:(numb_of_patients+4)]
greedy.ys <- greedy.results[1:(length(greedy.results)-4)]
greedy.power <- greedy.results[length(greedy.results)-3]
greedy.ENS <- mean(greedy.ys)
greedy.S10 <- greedy.results[length(greedy.results)]
sprintf("power is %f", greedy.power)

sim_data <- simulator_new(numb_of_sims, numb_of_patients, one_big_run, greedy, NA)
greedy.ys.1 <- sim_data[1:numb_of_patients]
greedy.props.1 <- sim_data[(numb_of_patients+5):(numb_of_sims+numb_of_patients+5)]
greedy.power.1 <- sim_data[numb_of_patients+5]
greedy.ENS.1 <- mean(greedy.ys)
greedy.S10.1 <- sim_data[numb_of_patients+4]


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper, NA)
thomper.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper.results <- sim_data[1:(numb_of_patients+4)]
thomper.ys <- thomper.results[1:(length(thomper.results)-4)]
thomper.power <- thomper.results[length(thomper.results)-3]
thomper.ENS <- mean(thomper.ys)
thomper.S10 <- thomper.results[length(thomper.results)]
sprintf("power is %f", thomper.power)

if (FALSE) {
sim_data <- simulator_new(numb_of_sims, numb_of_patients, one_big_run, thomper, NA)
thomper.ys.1 <- sim_data[1:numb_of_patients]
thomper.power.1 <- sim_data[numb_of_patients+1]
thomper.S10.1 <- sim_data[numb_of_patients+4]
thomper.ENS.1 <- mean(thomper.ys)
thomper.props.1 <- sim_data[(numb_of_patients+5):(numb_of_sims+numb_of_patients+5)]
}

if (FALSE) {
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw_foo, 1)
thomper_foo.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_foo.results <- sim_data[1:(numb_of_patients+4)]
thomper_foo.ys <- thomper_foo.results[1:(length(thomper_foo.results)-4)]
thomper_foo.power <- thomper_foo.results[length(thomper_foo.results)-3]
thomper_foo.ENS <- mean(thomper_foo.ys)
thomper_foo.S10 <- thomper_foo.results[length(thomper_foo.results)]
sprintf("power is %f", thomper_foo.power)
}
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, balanced, NA)
balanced.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
balanced.results <- sim_data[1:(numb_of_patients+4)]
balanced.ys <- balanced.results[1:(length(balanced.results)-4)]
balanced.power <- balanced.results[length(balanced.results)-3]
balanced.ENS <- mean(balanced.ys)
balanced.S10 <- balanced.results[length(balanced.results)]
sprintf("power is %f", balanced.power)

if (FALSE) {
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 0.5)
tw.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw.results <- sim_data[1:(numb_of_patients+4)]
tw.ys <- tw.results[1:(length(tw.results)-4)]
tw.power <- tw.results[length(tw.results)-3]
tw.ENS <- mean(tw.ys)
tw.S10 <- tw.results[length(tw.results)]
sprintf("power is %f", tw.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 1/9)
tw_one_ninth.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_one_ninth.results <- sim_data[1:(numb_of_patients+4)]
tw_one_ninth.ys <- tw_one_ninth.results[1:(length(tw_one_ninth.results)-4)]
tw_one_ninth.power <- tw_one_ninth.results[length(tw_one_ninth.results)-3]
tw_one_ninth.ENS <- mean(tw_one_ninth.ys)
tw_one_ninth.S10 <- tw_one_ninth.results[length(tw_one_ninth.results)]
sprintf("power is %f", tw_one_ninth.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 1/5)
tw_one_fifth.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_one_fifth.results <- sim_data[1:(numb_of_patients+4)]
tw_one_fifth.ys <- tw_one_fifth.results[1:(length(tw_one_fifth.results)-4)]
tw_one_fifth.power <- tw_one_fifth.results[length(tw_one_fifth.results)-3]
tw_one_fifth.ENS <- mean(tw_one_fifth.ys)
tw_one_fifth.S10 <- tw_one_fifth.results[length(tw_one_fifth.results)]
sprintf("power is %f", tw_one_fifth.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 1/3)
tw_one_third.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_one_third.results <- sim_data[1:(numb_of_patients+4)]
tw_one_third.ys <- tw_one_third.results[1:(length(tw_one_third.results)-4)]
tw_one_third.power <- tw_one_third.results[length(tw_one_third.results)-3]
tw_one_third.ENS <- mean(tw_one_third.ys)
tw_one_third.S10 <- tw_one_third.results[length(tw_one_third.results)]
sprintf("power is %f", tw_one_third.power)


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 2/3)
tw_two_thirds.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_two_thirds.results <- sim_data[1:(numb_of_patients+4)]
tw_two_thirds.ys <- tw_two_thirds.results[1:(length(tw_two_thirds.results)-4)]
tw_two_thirds.power <- tw_two_thirds.results[length(tw_two_thirds.results)-3]
tw_two_thirds.ENS <- mean(tw_two_thirds.ys)
tw_two_thirds.S10 <- tw_two_thirds.results[length(tw_two_thirds.results)]
sprintf("power is %f", tw_two_thirds.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 4/5)
tw_four_fifths.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_four_fifths.results <- sim_data[1:(numb_of_patients+4)]
tw_four_fifths.ys <- tw_four_fifths.results[1:(length(tw_four_fifths.results)-4)]
tw_four_fifths.power <- tw_four_fifths.results[length(tw_four_fifths.results)-3]
tw_four_fifths.ENS <- mean(tw_four_fifths.ys)
tw_four_fifths.S10 <- tw_four_fifths.results[length(tw_four_fifths.results)]
sprintf("power is %f", tw_four_fifths.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw, 8/9)
tw_eight_ninths.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_eight_ninths.results <- sim_data[1:(numb_of_patients+4)]
tw_eight_ninths.ys <- tw_eight_ninths.results[1:(length(tw_eight_ninths.results)-4)]
tw_eight_ninths.power <- tw_eight_ninths.results[length(tw_eight_ninths.results)-3]
tw_eight_ninths.ENS <- mean(tw_eight_ninths.ys)
tw_eight_ninths.S10 <- tw_eight_ninths.results[length(tw_eight_ninths.results)]
sprintf("power is %f", tw_eight_ninths.power)
}


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 1/9)
thomper_interperlated_one_ninth.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_one_ninth.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_one_ninth.ys <- thomper_interperlated_one_ninth.results[1:(length(thomper_interperlated_one_ninth.results)-4)]
thomper_interperlated_one_ninth.power <- thomper_interperlated_one_ninth.results[length(thomper_interperlated_one_ninth.results)-3]
thomper_interperlated_one_ninth.ENS <- mean(thomper_interperlated_one_ninth.ys)
thomper_interperlated_one_ninth.S10 <- thomper_interperlated_one_ninth.results[length(thomper_interperlated_one_ninth.results)]
sprintf("power is %f", thomper_interperlated_one_ninth.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 1/5)
thomper_interperlated_one_fifth.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_one_fifth.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_one_fifth.ys <- thomper_interperlated_one_fifth.results[1:(length(thomper_interperlated_one_fifth.results)-4)]
thomper_interperlated_one_fifth.power <- thomper_interperlated_one_fifth.results[length(thomper_interperlated_one_fifth.results)-3]
thomper_interperlated_one_fifth.ENS <- mean(thomper_interperlated_one_fifth.ys)
thomper_interperlated_one_fifth.S10 <- thomper_interperlated_one_fifth.results[length(thomper_interperlated_one_fifth.results)]
sprintf("power is %f", thomper_interperlated_one_fifth.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 1/3)
thomper_interperlated_one_third.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_one_third.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_one_third.ys <- thomper_interperlated_one_third.results[1:(length(thomper_interperlated_one_third.results)-4)]
thomper_interperlated_one_third.power <- thomper_interperlated_one_third.results[length(thomper_interperlated_one_third.results)-3]
thomper_interperlated_one_third.ENS <- mean(thomper_interperlated_one_third.ys)
thomper_interperlated_one_third.S10 <- thomper_interperlated_one_third.results[length(thomper_interperlated_one_third.results)]
sprintf("power is %f", thomper_interperlated_one_third.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 1/2)
thomper_interperlated_half.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_half.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_half.ys <- thomper_interperlated_half.results[1:(length(thomper_interperlated_half.results)-4)]
thomper_interperlated_half.power <- thomper_interperlated_half.results[length(thomper_interperlated_half.results)-3]
thomper_interperlated_half.ENS <- mean(thomper_interperlated_half.ys)
thomper_interperlated_half.S10 <- thomper_interperlated_half.results[length(thomper_interperlated_half.results)]
sprintf("power is %f", thomper_interperlated_half.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 2/3)
thomper_interperlated_two_thirds.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_two_thirds.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_two_thirds.ys <- thomper_interperlated_two_thirds.results[1:(length(thomper_interperlated_two_thirds.results)-4)]
thomper_interperlated_two_thirds.power <- thomper_interperlated_two_thirds.results[length(thomper_interperlated_two_thirds.results)-3]
thomper_interperlated_two_thirds.ENS <- mean(thomper_interperlated_two_thirds.ys)
thomper_interperlated_two_thirds.S10 <- thomper_interperlated_two_thirds.results[length(thomper_interperlated_two_thirds.results)]
sprintf("power is %f", thomper_interperlated_two_thirds.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 4/5)
thomper_interperlated_four_fifths.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_four_fifths.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_four_fifths.ys <- thomper_interperlated_four_fifths.results[1:(length(thomper_interperlated_four_fifths.results)-4)]
thomper_interperlated_four_fifths.power <- thomper_interperlated_four_fifths.results[length(thomper_interperlated_four_fifths.results)-3]
thomper_interperlated_four_fifths.ENS <- mean(thomper_interperlated_four_fifths.ys)
thomper_interperlated_four_fifths.S10 <- thomper_interperlated_four_fifths.results[length(thomper_interperlated_four_fifths.results)]
sprintf("power is %f", thomper_interperlated_four_fifths.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_interperlated, 8/9)
thomper_interperlated_eight_ninths.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_interperlated_eight_ninths.results <- sim_data[1:(numb_of_patients+4)]
thomper_interperlated_eight_ninths.ys <- thomper_interperlated_eight_ninths.results[1:(length(thomper_interperlated_eight_ninths.results)-4)]
thomper_interperlated_eight_ninths.power <- thomper_interperlated_eight_ninths.results[length(thomper_interperlated_eight_ninths.results)-3]
thomper_interperlated_eight_ninths.ENS <- mean(thomper_interperlated_eight_ninths.ys)
thomper_interperlated_eight_ninths.S10 <- thomper_interperlated_eight_ninths.results[length(thomper_interperlated_eight_ninths.results)]
sprintf("power is %f", thomper_interperlated_eight_ninths.power)


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,5,0))
thomper_add_5.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_add_5.results <- sim_data[1:(numb_of_patients+4)]
thomper_add_5.ys <- thomper_add_5.results[1:(length(thomper_add_5.results)-4)]
thomper_add_5.power <- thomper_add_5.results[length(thomper_add_5.results)-3]
thomper_add_5.ENS <- mean(thomper_add_5.ys)
thomper_add_5.S10 <- thomper_add_5.results[length(thomper_add_5.results)]
sprintf("power is %f", thomper_add_5.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,25,0))
thomper_add_25.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_add_25.results <- sim_data[1:(numb_of_patients+4)]
thomper_add_25.ys <- thomper_add_25.results[1:(length(thomper_add_25.results)-4)]
thomper_add_25.power <- thomper_add_25.results[length(thomper_add_25.results)-3]
thomper_add_25.ENS <- mean(thomper_add_25.ys)
thomper_add_25.S10 <- thomper_add_25.results[length(thomper_add_25.results)]
sprintf("power is %f", thomper_add_25.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,100,0))
thomper_add_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_add_100.results <- sim_data[1:(numb_of_patients+4)]
thomper_add_100.ys <- thomper_add_100.results[1:(length(thomper_add_100.results)-4)]
thomper_add_100.power <- thomper_add_100.results[length(thomper_add_100.results)-3]
thomper_add_100.ENS <- mean(thomper_add_100.ys)
thomper_add_100.S10 <- thomper_add_100.results[length(thomper_add_100.results)]
sprintf("power is %f", thomper_add_100.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(0.5,0,0))
thomper_divide_2.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_divide_2.results <- sim_data[1:(numb_of_patients+4)]
thomper_divide_2.ys <- thomper_divide_2.results[1:(length(thomper_divide_2.results)-4)]
thomper_divide_2.power <- thomper_divide_2.results[length(thomper_divide_2.results)-3]
thomper_divide_2.ENS <- mean(thomper_divide_2.ys)
thomper_divide_2.S10 <- thomper_divide_2.results[length(thomper_divide_2.results)]
sprintf("power is %f", thomper_divide_2.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(0.2,0,0))
thomper_divide_5.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_divide_5.results <- sim_data[1:(numb_of_patients+4)]
thomper_divide_5.ys <- thomper_divide_5.results[1:(length(thomper_divide_5.results)-4)]
thomper_divide_5.power <- thomper_divide_5.results[length(thomper_divide_5.results)-3]
thomper_divide_5.ENS <- mean(thomper_divide_5.ys)
thomper_divide_5.S10 <- thomper_divide_5.results[length(thomper_divide_5.results)]
sprintf("power is %f", thomper_divide_5.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1/25,0,0))
thomper_divide_25.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_divide_25.results <- sim_data[1:(numb_of_patients+4)]
thomper_divide_25.ys <- thomper_divide_25.results[1:(length(thomper_divide_25.results)-4)]
thomper_divide_25.power <- thomper_divide_25.results[length(thomper_divide_25.results)-3]
thomper_divide_25.ENS <- mean(thomper_divide_25.ys)
thomper_divide_25.S10 <- thomper_divide_25.results[length(thomper_divide_25.results)]
sprintf("power is %f", thomper_divide_25.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1/100,0,0))
thomper_divide_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_divide_100.results <- sim_data[1:(numb_of_patients+4)]
thomper_divide_100.ys <- thomper_divide_100.results[1:(length(thomper_divide_100.results)-4)]
thomper_divide_100.power <- thomper_divide_100.results[length(thomper_divide_100.results)-3]
thomper_divide_100.ENS <- mean(thomper_divide_100.ys)
thomper_divide_100.S10 <- thomper_divide_100.results[length(thomper_divide_100.results)]
sprintf("power is %f", thomper_divide_100.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(2,0,0))
thomper_multiply_2.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_multiply_2.results <- sim_data[1:(numb_of_patients+4)]
thomper_multiply_2.ys <- thomper_multiply_2.results[1:(length(thomper_multiply_2.results)-4)]
thomper_multiply_2.power <- thomper_multiply_2.results[length(thomper_multiply_2.results)-3]
thomper_multiply_2.ENS <- mean(thomper_multiply_2.ys)
thomper_multiply_2.S10 <- thomper_multiply_2.results[length(thomper_multiply_2.results)]
sprintf("power is %f", thomper_multiply_2.power)


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_root, 0.5)
thomper_sqrt.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_sqrt.results <- sim_data[1:(numb_of_patients+4)]
thomper_sqrt.ys <- thomper_sqrt.results[1:(length(thomper_sqrt.results)-4)]
thomper_sqrt.power <- thomper_sqrt.results[length(thomper_sqrt.results)-3]
thomper_sqrt.ENS <- mean(thomper_sqrt.ys)
thomper_sqrt.S10 <- thomper_sqrt.results[length(thomper_sqrt.results)]
sprintf("power is %f", thomper_sqrt.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_root, 1/3)
thomper_cbrt.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_cbrt.results <- sim_data[1:(numb_of_patients+4)]
thomper_cbrt.ys <- thomper_cbrt.results[1:(length(thomper_cbrt.results)-4)]
thomper_cbrt.power <- thomper_cbrt.results[length(thomper_cbrt.results)-3]
thomper_cbrt.ENS <- mean(thomper_cbrt.ys)
thomper_cbrt.S10 <- thomper_cbrt.results[length(thomper_cbrt.results)]
sprintf("power is %f", thomper_cbrt.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_log, exp(1))
thomper_log_e.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_log_e.results <- sim_data[1:(numb_of_patients+4)]
thomper_log_e.ys <- thomper_log_e.results[1:(length(thomper_log_e.results)-4)]
thomper_log_e.power <- thomper_log_e.results[length(thomper_log_e.results)-3]
thomper_log_e.ENS <- mean(thomper_log_e.ys)
thomper_log_e.S10 <- thomper_log_e.results[length(thomper_log_e.results)]
sprintf("power is %f", thomper_log_e.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_log, 10)
thomper_log_10.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_log_10.results <- sim_data[1:(numb_of_patients+4)]
thomper_log_10.ys <- thomper_log_10.results[1:(length(thomper_log_10.results)-4)]
thomper_log_10.power <- thomper_log_10.results[length(thomper_log_10.results)-3]
thomper_log_10.ENS <- mean(thomper_log_10.ys)
thomper_log_10.S10 <- thomper_log_10.results[length(thomper_log_10.results)]
sprintf("power is %f", thomper_log_10.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,-5,1))
thomper_subtract_5.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_subtract_5.results <- sim_data[1:(numb_of_patients+4)]
thomper_subtract_5.ys <- thomper_subtract_5.results[1:(length(thomper_subtract_5.results)-4)]
thomper_subtract_5.power <- thomper_subtract_5.results[length(thomper_subtract_5.results)-3]
thomper_subtract_5.ENS <- mean(thomper_subtract_5.ys)
thomper_subtract_5.S10 <- thomper_subtract_5.results[length(thomper_subtract_5.results)]
sprintf("power is %f", thomper_subtract_5.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,-25,1))
thomper_subtract_25.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_subtract_25.results <- sim_data[1:(numb_of_patients+4)]
thomper_subtract_25.ys <- thomper_subtract_25.results[1:(length(thomper_subtract_25.results)-4)]
thomper_subtract_25.power <- thomper_subtract_25.results[length(thomper_subtract_25.results)-3]
thomper_subtract_25.ENS <- mean(thomper_subtract_25.ys)
thomper_subtract_25.S10 <- thomper_subtract_25.results[length(thomper_subtract_25.results)]
sprintf("power is %f", thomper_subtract_25.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,-50,1))
thomper_subtract_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_subtract_100.results <- sim_data[1:(numb_of_patients+4)]
thomper_subtract_100.ys <- thomper_subtract_100.results[1:(length(thomper_subtract_100.results)-4)]
thomper_subtract_100.power <- thomper_subtract_100.results[length(thomper_subtract_100.results)-3]
thomper_subtract_100.ENS <- mean(thomper_subtract_100.ys)
thomper_subtract_100.S10 <- thomper_subtract_100.results[length(thomper_subtract_100.results)]
sprintf("power is %f", thomper_subtract_100.power)


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,0,5))
thomper_min_5.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_min_5.results <- sim_data[1:(numb_of_patients+4)]
thomper_min_5.ys <- thomper_min_5.results[1:(length(thomper_min_5.results)-4)]
thomper_min_5.power <- thomper_min_5.results[length(thomper_min_5.results)-3]
thomper_min_5.ENS <- mean(thomper_min_5.ys)
thomper_min_5.S10 <- thomper_min_5.results[length(thomper_min_5.results)]
sprintf("power is %f", thomper_min_5.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,0,25))
thomper_min_25.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_min_25.results <- sim_data[1:(numb_of_patients+4)]
thomper_min_25.ys <- thomper_min_25.results[1:(length(thomper_min_25.results)-4)]
thomper_min_25.power <- thomper_min_25.results[length(thomper_min_25.results)-3]
thomper_min_25.ENS <- mean(thomper_min_25.ys)
thomper_min_25.S10 <- thomper_min_25.results[length(thomper_min_25.results)]
sprintf("power is %f", thomper_min_25.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_trans, c(1,0,100))
thomper_min_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_min_100.results <- sim_data[1:(numb_of_patients+4)]
thomper_min_100.ys <- thomper_min_100.results[1:(length(thomper_min_100.results)-4)]
thomper_min_100.power <- thomper_min_100.results[length(thomper_min_100.results)-3]
thomper_min_100.ENS <- mean(thomper_min_100.ys)
thomper_min_100.S10 <- thomper_min_100.results[length(thomper_min_100.results)]
sprintf("power is %f", thomper_min_100.power)

if (FALSE) {
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw_foo, 0.5)
tw_foo.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_foo.results <- sim_data[1:(numb_of_patients+4)]
tw_foo.ys <- tw_foo.results[1:(length(tw_foo.results)-4)]
tw_foo.power <- tw_foo.results[length(tw_foo.results)-3]
tw_foo.ENS <- mean(tw_foo.ys)
tw_foo.S10 <- tw_foo.results[length(tw_foo.results)]
sprintf("power is %f", tw_foo.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, tw_falling, NA)
tw_falling.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
tw_falling.results <- sim_data[1:(numb_of_patients+4)]
tw_falling.ys <- tw_falling.results[1:(length(tw_falling.results)-4)]
tw_falling.power <- tw_falling.results[length(tw_falling.results)-3]
tw_falling.ENS <- mean(tw_falling.ys)
tw_falling.S10 <- tw_falling.results[length(tw_falling.results)]
sprintf("power is %f", tw_falling.power)
}
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, Neyman_with_updating, NA)
Neyman_with_updating.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
Neyman_with_updating.results <- sim_data[1:(numb_of_patients+4)]
Neyman_with_updating.ys <- Neyman_with_updating.results[1:(length(Neyman_with_updating.results)-4)]
Neyman_with_updating.power <- Neyman_with_updating.results[length(Neyman_with_updating.results)-3]
Neyman_with_updating.ENS <- mean(Neyman_with_updating.ys)
Neyman_with_updating.S10 <- Neyman_with_updating.results[length(Neyman_with_updating.results)]
sprintf("power is %f", Neyman_with_updating.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, rosey, NA)
rosey.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
rosey.results <- sim_data[1:(numb_of_patients+4)]
rosey.ys <- rosey.results[1:(length(rosey.results)-4)]
rosey.power <- rosey.results[length(rosey.results)-3]
rosey.ENS <- mean(rosey.ys)
rosey.S10 <- rosey.results[length(rosey.results)]
sprintf("power is %f", rosey.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, DTL, 0)
DTL_2.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
DTL_2.results <- sim_data[1:(numb_of_patients+4)]
DTL_2.ys <- DTL_2.results[1:(length(DTL_2.results)-4)]
DTL_2.power <- DTL_2.results[length(DTL_2.results)-3]
DTL_2.ENS <- mean(DTL_2.ys)
DTL_2.S10 <- DTL_2.results[length(DTL_2.results)]
sprintf("power is %f", DTL_2.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, DTL, 2)
DTL_4.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
DTL_4.results <- sim_data[1:(numb_of_patients+4)]
DTL_4.ys <- DTL_4.results[1:(length(DTL_4.results)-4)]
DTL_4.power <- DTL_4.results[length(DTL_4.results)-3]
DTL_4.ENS <- mean(DTL_4.ys)
DTL_4.S10 <- DTL_4.results[length(DTL_4.results)]
sprintf("power is %f", DTL_4.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, DTL, 8)
DTL_10.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
DTL_10.results <- sim_data[1:(numb_of_patients+4)]
DTL_10.ys <- DTL_10.results[1:(length(DTL_10.results)-4)]
DTL_10.power <- DTL_10.results[length(DTL_10.results)-3]
DTL_10.ENS <- mean(DTL_10.ys)
DTL_10.S10 <- DTL_10.results[length(DTL_10.results)]
sprintf("power is %f", DTL_10.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, DBCD, 0)
DBCD_0.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
DBCD_0.results <- sim_data[1:(numb_of_patients+4)]
DBCD_0.ys <- DBCD_0.results[1:(length(DBCD_0.results)-4)]
DBCD_0.power <- DBCD_0.results[length(DBCD_0.results)-3]
DBCD_0.ENS <- mean(DBCD_0.ys)
DBCD_0.S10 <- DBCD_0.results[length(DBCD_0.results)]
sprintf("power is %f", DBCD_0.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, DBCD, 2)
DBCD_2.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
DBCD_2.results <- sim_data[1:(numb_of_patients+4)]
DBCD_2.ys <- DBCD_2.results[1:(length(DBCD_2.results)-4)]
DBCD_2.power <- DBCD_2.results[length(DBCD_2.results)-3]
DBCD_2.ENS <- mean(DBCD_2.ys)
DBCD_2.S10 <- DBCD_2.results[length(DBCD_2.results)]
sprintf("power is %f", DBCD_2.power)

if (FALSE) {
numb_of_sims <- 100
sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, FLGI, 2)
FLGI_2.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
FLGI_2.results <- sim_data[1:(numb_of_patients+4)]
FLGI_2.ys <- FLGI_2.results[1:(length(FLGI_2.results)-4)]
FLGI_2.power <- FLGI_2.results[length(FLGI_2.results)-3]
FLGI_2.ENS <- mean(FLGI_2.ys)
FLGI_2.S10 <- FLGI_2.results[length(FLGI_2.results)]
sprintf("power is %f", FLGI_2.power)


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, FLGI, 4)
FLGI_4.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
FLGI_4.results <- sim_data[1:(numb_of_patients+4)]
FLGI_4.ys <- FLGI_4.results[1:(length(FLGI_4.results)-4)]
FLGI_4.power <- FLGI_4.results[length(FLGI_4.results)-3]
FLGI_4.ENS <- mean(FLGI_4.ys)
FLGI_4.S10 <- FLGI_4.results[length(FLGI_4.results)]
sprintf("power is %f", FLGI_4.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, FLGI, 6)
FLGI_6.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
FLGI_6.results <- sim_data[1:(numb_of_patients+4)]
FLGI_6.ys <- FLGI_6.results[1:(length(FLGI_6.results)-4)]
FLGI_6.power <- FLGI_6.results[length(FLGI_6.results)-3]
FLGI_6.ENS <- mean(FLGI_6.ys)
FLGI_6.S10 <- FLGI_6.results[length(FLGI_6.results)]
sprintf("power is %f", FLGI_6.power)
}


sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, ebm, 100)
ebm_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
ebm_100.results <- sim_data[1:(numb_of_patients+4)]
ebm_100.ys <- ebm_100.results[1:(length(ebm_100.results)-4)]
ebm_100.power <- ebm_100.results[length(ebm_100.results)-3]
ebm_100.ENS <- mean(ebm_100.ys)
ebm_100.S10 <- ebm_100.results[length(ebm_100.results)]
sprintf("power is %f", ebm_100.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, ebm, 1000)
ebm_1000.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
ebm_1000.results <- sim_data[1:(numb_of_patients+4)]
ebm_1000.ys <- ebm_1000.results[1:(length(ebm_1000.results)-4)]
ebm_1000.power <- ebm_1000.results[length(ebm_1000.results)-3]
ebm_1000.ENS <- mean(ebm_1000.ys)
ebm_1000.S10 <- ebm_1000.results[length(ebm_1000.results)]
sprintf("power is %f", ebm_1000.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_falling_min, c(100,1))
thomper_falling_min_100.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_falling_min_100.results <- sim_data[1:(numb_of_patients+4)]
thomper_falling_min_100.ys <- thomper_falling_min_100.results[1:(length(thomper_falling_min_100.results)-4)]
thomper_falling_min_100.power <- thomper_falling_min_100.results[length(thomper_falling_min_100.results)-3]
thomper_falling_min_100.ENS <- mean(thomper_falling_min_100.ys)
thomper_falling_min_100.S10 <- thomper_falling_min_100.results[length(thomper_falling_min_100.results)]
sprintf("power is %f", thomper_falling_min_100.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_falling_min, c(1000,1))
thomper_falling_min_1000.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_falling_min_1000.results <- sim_data[1:(numb_of_patients+4)]
thomper_falling_min_1000.ys <- thomper_falling_min_1000.results[1:(length(thomper_falling_min_1000.results)-4)]
thomper_falling_min_1000.power <- thomper_falling_min_1000.results[length(thomper_falling_min_1000.results)-3]
thomper_falling_min_1000.ENS <- mean(thomper_falling_min_1000.ys)
thomper_falling_min_1000.S10 <- thomper_falling_min_1000.results[length(thomper_falling_min_1000.results)]
sprintf("power is %f", thomper_falling_min_1000.power)

sim_data <- simulator(numb_of_sims, numb_of_patients, one_run, thomper_falling_min, c(10000,1))
thomper_falling_min_10000.props <- sim_data[(numb_of_patients+6):(length(sim_data))]
thomper_falling_min_10000.results <- sim_data[1:(numb_of_patients+4)]
thomper_falling_min_10000.ys <- thomper_falling_min_10000.results[1:(length(thomper_falling_min_10000.results)-4)]
thomper_falling_min_10000.power <- thomper_falling_min_10000.results[length(thomper_falling_min_10000.results)-3]
thomper_falling_min_10000.ENS <- mean(thomper_falling_min_10000.ys)
thomper_falling_min_10000.S10 <- thomper_falling_min_10000.results[length(thomper_falling_min_10000.results)]
sprintf("power is %f", thomper_falling_min_10000.power)


## plot hist of proportions

# greedy.hist <- hist(greedy.props,plot=FALSE,breaks=breaks)
# thomper.hist <- hist(thomper.props,plot=FALSE,breaks=breaks)
# balanced.hist <- hist(balanced.props,plot=FALSE,breaks=breaks)
# tw.hist <- hist(tw.props,plot=FALSE,breaks=breaks)
# DTL.hist <- hist(DTL.props,plot=FALSE,breaks=breaks)

# #plot(hist1, col = rgb(1,0,0,0.4),freq = TRUE)
# plot(hist2,add=FALSE,col = rgb(0,0,1,0.4),freq = FALSE)
# # plot(hist3,add=FALSE,col = rgb(0,1,1,0.4),freq = FALSE)
# plot(hist4,add=TRUE,col = rgb(0,1,0,0.4),freq = FALSE)
# plot(DTL.hist,add=TRUE,col = rgb(1,0,0,0.4),freq = FALSE)
# plot(hist(thomper.props),add=TRUE,col = rgb(0,0,1,0.4))




## plot sim progressions

x <- 1:(numb_of_patients)
colours <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
bandits <- c("Greedy","TS","TS/25","DTL(2)","TS Falling Min 10,000","TS/100","Rosey","TS Falling Min 1000","EBM n=1000")
width <- 1

plot(greedy.ys ~ x, type='l',xlim=c(0,100000),ylim=c(0,1), main="Comparison of Systems", sub="100,000 Patients", xlab="Patient Index", ylab="Probability of Getting the Best Treatment", col=colours[1], lwd=width)

#lines(thomper.ys ~ x, col="red")
fit <- glm(formula = thomper.ys ~ x * log(x) * log(x^2), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[2], lwd=width)

#lines(thomper_divide_25.ys ~ x, col="blue")
fit <- glm(formula = thomper_divide_25.ys ~ x * log(x) * log(x^2), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[3], lwd=width)

#lines(DTL_2.ys ~ x, col="black")
fit <- glm(formula = DTL_2.ys ~ x * log(x) * log(x^2), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[4], lwd=width)

#lines(DBCD_0.ys ~ x, col="green")
segments(x0=0,y0=mean(thomper_falling_min_10000.ys[1:5000]),x1=5000,y1=mean(thomper_falling_min_10000.ys[1:5000]))
thomper_falling_min_10000.ys <- thomper_falling_min_10000.ys[5001:numb_of_patients]
fit <- glm(formula = thomper_falling_min_10000.ys ~ x * x^2 * x^3, family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[5], lwd=width)
lines(thomper_falling_min_10000.ys ~ x)
abline(v=5000)

#lines(thomper_divide_100.ys ~ x, col="purple")
fit <- glm(formula = thomper_divide_100.ys ~ x * log(x) * log(x^2), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[6], lwd=width)

#lines(rosey.ys ~ x, col="brown")
fit <- glm(formula = rosey.ys ~ x * log(x) * log(x^2), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[7], lwd=width)

#lines(thomper_min_100.ys ~ x, col="orange")
fit <- glm(formula = thomper_falling_min_1000.ys ~ x * log(x), family = binomial)
newdat <- data.frame(ys=x)
newdat$ys <- predict(fit, newdata=newdat, type="response")
lines(ys ~ x, newdat, col=colours[8], lwd=width)


#lines(thomper_add_100.ys ~ x, col="dark-grey")
fit <- glm(formula = ebm_1000.ys ~ x * log(x) * log(x^2), family = binomial)
segments(x0=0, y0=mean(ebm_1000.ys[1:1000]), x1=1000, y1=mean(ebm_1000.ys[1:1000]),col=colours[9],lwd=width)
segments(x0=1001, y0=mean(ebm_1000.ys[1001:numb_of_patients]), x1=100000, y1=mean(ebm_1000.ys[1001:numb_of_patients]),col=colours[9],lwd=width)
segments(x0=1000, y0=mean(ebm_1000.ys[1:1000]), x1=1001, y1=mean(ebm_1000.ys[1001:numb_of_patients]),col=colours[9],lwd=width)


legend(x="bottomright", legend=bandits, col=colours, pch=15)


# lines(thomper_foo.ys ~ x, col="purple")
# lines(balanced.ys ~ x, col="blue")
# lines(tw.ys ~ x, col="yellow")
# lines(tw_foo.ys ~ x, col="black")


## plot tradeoff between power and ENS
x <- c(greedy.power,thomper.power,balanced.power,tw.power,Neyman_with_updating.power,rosey.power,tw_one_ninth.power,tw_one_fifth.power,tw_one_third.power,tw_two_thirds.power,tw_four_fifths.power,tw_eight_ninths.power,DTL_2.power,DTL_4.power,DTL_10.power,DBCD_0.power,DBCD_2.power,thomper_interperlated_one_ninth.power,thomper_interperlated_one_fifth.power,thomper_interperlated_one_third.power,thomper_interperlated_half.power,thomper_interperlated_two_thirds.power,thomper_interperlated_four_fifths.power,thomper_interperlated_eight_ninths.power,thomper_add_5.power,thomper_add_25.power,thomper_add_100.power,thomper_divide_2.power,thomper_divide_5.power,thomper_divide_25.power,thomper_divide_100.power,thomper_multiply_2.power,thomper_sqrt.power,thomper_cbrt.power,thomper_log_e.power,thomper_log_10.power,thomper_subtract_5.power,thomper_subtract_25.power,thomper_subtract_100.power,thomper_min_5.power,thomper_min_25.power,thomper_min_100.power)#,FLGI_2.power)#,FLGI_4.power,FLGI_6.power)
y <- c(greedy.ENS,  thomper.ENS,  balanced.ENS,  tw.ENS,  Neyman_with_updating.ENS,  rosey.ENS,  tw_one_ninth.ENS,  tw_one_fifth.ENS,  tw_one_third.ENS,  tw_two_thirds.ENS,  tw_four_fifths.ENS,  tw_eight_ninths.ENS,  DTL_2.ENS,  DTL_4.ENS,  DTL_10.ENS,  DBCD_0.ENS,  DBCD_2.ENS,  thomper_interperlated_one_ninth.ENS,  thomper_interperlated_one_fifth.ENS,  thomper_interperlated_one_third.ENS,  thomper_interperlated_half.ENS,  thomper_interperlated_two_thirds.ENS,  thomper_interperlated_four_fifths.ENS,  thomper_interperlated_eight_ninths.ENS,  thomper_add_5.ENS,  thomper_add_25.ENS,  thomper_add_100.ENS,  thomper_divide_2.ENS,  thomper_divide_5.ENS,  thomper_divide_25.ENS,  thomper_divide_100.ENS,  thomper_multiply_2.ENS,  thomper_sqrt.ENS,  thomper_cbrt.ENS,  thomper_log_e.ENS,  thomper_log_10.ENS,  thomper_subtract_5.ENS,  thomper_subtract_25.ENS,  thomper_subtract_100.ENS,  thomper_min_5.ENS,  thomper_min_25.ENS,  thomper_min_100.ENS  )#,FLGI_2.ENS)#,FLGI_4.ENS,FLGI_6.ENS)
labels <- c("Greedy","Thompson Sampling","Equal Randomisation","TW(1/2)","Neyman",   "Rosey",    "TW(1/9)",         "TW(1/5)",         "TW(1/3)",         "TW(2/3)",          "TW(4/5)",           "TW(8/9)",            "DTL(2)",   "DTL(4)",   "DTL(10)",   "DBCD(0)",   "DBCD(2)",   "TI(1/9)",                            "TI(1/5)",                            "TI(1/3)",                            "TI(1/2)",                       "TI(2/3)",                             "TI(4/5)",                              "TI(8/9)",                               "TS+5",             "TS+25",             "TS+100",             "TS/2",                "TS/5",                "TS/25",                "TS/100",                "TS*2",                  "TS_sqrt",         "TS_cbrt",         "TS_ln",            "TS_log_10",         "TS-5",                  "TS-25",                  "TS-100",                  "TS_min_5",         "TS_min_25",         "TS_min_100")#,"FLGI(2)")#,"FLGI(4)","FLGI(6)")

x_main <- c(greedy.power,thomper.power,balanced.power,tw.power, tw_one_fifth.power, DTL_10.power)#,FLGI_2.power)
y_main <- c(greedy.ENS,thomper.ENS,balanced.ENS,tw.ENS, tw_one_fifth.ENS, DTL_10.ENS)#,FLGI_2.ENS)
labels_main <- c("Greedy","Thompson Sampling","Equal Randomisation","TW(1/2)", "TW(1/5)", "DTL(10)")#,"FLGI(2)")

x_fast <- c(greedy.power,thomper.power,balanced.power,Neyman_with_updating.power,rosey.power,DTL_2.power,DTL_4.power,DTL_10.power,DBCD_0.power,DBCD_2.power,thomper_interperlated_one_ninth.power,thomper_interperlated_one_fifth.power,thomper_interperlated_one_third.power,thomper_interperlated_half.power,thomper_interperlated_two_thirds.power,thomper_interperlated_four_fifths.power,thomper_interperlated_eight_ninths.power,thomper_add_5.power,thomper_add_25.power,thomper_add_100.power,thomper_divide_2.power,thomper_divide_5.power,thomper_divide_25.power,thomper_divide_100.power,thomper_multiply_2.power,thomper_sqrt.power,thomper_cbrt.power,thomper_log_e.power,thomper_log_10.power,thomper_subtract_5.power,thomper_subtract_25.power,thomper_subtract_100.power,thomper_min_5.power,thomper_min_25.power,thomper_min_100.power)#,FLGI_2.power)#,FLGI_4.power,FLGI_6.power)
y_fast <- c(greedy.ENS,  thomper.ENS,  balanced.ENS,  Neyman_with_updating.ENS,  rosey.ENS,  DTL_2.ENS,  DTL_4.ENS,  DTL_10.ENS,  DBCD_0.ENS,  DBCD_2.ENS,  thomper_interperlated_one_ninth.ENS,  thomper_interperlated_one_fifth.ENS,  thomper_interperlated_one_third.ENS,  thomper_interperlated_half.ENS,  thomper_interperlated_two_thirds.ENS,  thomper_interperlated_four_fifths.ENS,  thomper_interperlated_eight_ninths.ENS,  thomper_add_5.ENS,  thomper_add_25.ENS,  thomper_add_100.ENS,  thomper_divide_2.ENS,  thomper_divide_5.ENS,  thomper_divide_25.ENS,  thomper_divide_100.ENS,  thomper_multiply_2.ENS,  thomper_sqrt.ENS,  thomper_cbrt.ENS,  thomper_log_e.ENS,  thomper_log_10.ENS,  thomper_subtract_5.ENS,  thomper_subtract_25.ENS,  thomper_subtract_100.ENS,  thomper_min_5.ENS,  thomper_min_25.ENS,  thomper_min_100.ENS  )#,FLGI_2.ENS)#,FLGI_4.ENS,FLGI_6.ENS)
labels_fast <- c("Greedy","Thompson Sampling","Equal Randomisation","Neyman",    "Rosey",    "DTL(2)",   "DTL(4)",   "DTL(10)",   "DBCD(0)",   "DBCD(2)",   "TI(1/9)",                            "TI(1/5)",                            "TI(1/3)",                            "TI(1/2)",                       "TI(2/3)",                             "TI(4/5)",                              "TI(8/9)",                               "TS+5",             "TS+25",             "TS+100",             "TS/2",                "TS/5",                "TS/25",                "TS/100",                "TS*2",                  "TS_sqrt",         "TS_cbrt",         "TS_ln",            "TS_log_10",         "TS-5",                  "TS-25",                  "TS-100",                  "TS_min_5",         "TS_min_25",         "TS_min_100")#,"FLGI(2)")#,"FLGI(4)","FLGI(6)")

x <- x_fast
y <- y_fast
labels <- labels_fast

plot(x_fast,y_fast,xlim=c(0.95,1),ylim=c(0.98,1),main="Expected Number of Successes vs Statistical Power",sub="100,000 Patients",xlab="Power",ylab="ENS")
#text(x+0.02,y+0.02,labels=labels,cex=0.5)
text(x_fast,x_fast,labels=labels_fast,cex=0.8)

# being on the convex hull isn't exactly what it means to be on the frontier.
# really the frontier is points where no other point is 'up and to the right' of it.
df <- data.frame(X=c(0,0,max(x_fast),x_fast), Y=c(0,max(y_fast),0,y_fast))
hull <- chull(df)
index_of_one <- match(1,hull)
hull <- append(hull[(index_of_one+1):length(hull)], hull[1:(index_of_one-1)])
hull <- df[hull,]
#hull <- rbind(df[hull,],df[hull[1],])
lines(hull)

tws <- data.frame(X=c(balanced.power,tw_one_ninth.power,tw_one_fifth.power,tw_one_third.power,tw.power,tw_two_thirds.power,tw_four_fifths.power,tw_eight_ninths.power,thomper.power), Y=c(balanced.ENS,tw_one_ninth.ENS,tw_one_fifth.ENS,tw_one_third.ENS,tw.ENS,tw_two_thirds.ENS,tw_four_fifths.ENS,tw_eight_ninths.ENS,thomper.ENS))
lines(tws, lty=1)

tis <- data.frame(X=c(balanced.power,thomper_interperlated_one_ninth.power,thomper_interperlated_one_fifth.power,thomper_interperlated_one_third.power,thomper_interperlated_half.power,thomper_interperlated_two_thirds.power,thomper_interperlated_four_fifths.power,thomper_interperlated_eight_ninths.power,thomper.power), Y=c(balanced.ENS,thomper_interperlated_one_ninth.ENS,thomper_interperlated_one_fifth.ENS,thomper_interperlated_one_third.ENS,thomper_interperlated_half.ENS,thomper_interperlated_two_thirds.ENS,thomper_interperlated_four_fifths.ENS,thomper_interperlated_eight_ninths.ENS,thomper.ENS))
lines(tis, lty=1)

dtls <- data.frame(X=c(DTL_2.power,DTL_4.power,DTL_10.power),Y=c(DTL_2.ENS,DTL_4.ENS,DTL_10.ENS))
lines(dtls, lty=2)

#flgis <- data.frame(X=c(FLGI_2.power,FLGI_4.power,FLGI_6.power),Y=c(FLGI_2.ENS,FLGI_4.ENS,FLGI_6.ENS))
#lines(flgis, lty=3)

# ## plot tradeoff between s10 and ENS

# plot(c(S101,S102,S104,S105,S106,S107,S108),c(ENS1,ENS2,ENS4,ENS5,ENS6,ENS7,ENS8),xlim=c(0,1),ylim=c(0,1),main="Expected Number of Successes vs s_0.1",xlab="s_0.1",ylab="ENS")
# text(c(S101,S102,S104,S105,S106,S107,S108)+0.07,c(ENS1,ENS2,ENS4,ENS5,ENS6,ENS7,ENS8)-0.00,labels=c("Greedy","Thomper","Balanced","TW(0.5)","TW(i/2n)","Neyman","Rosey"))


