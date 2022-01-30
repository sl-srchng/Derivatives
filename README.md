# Derivatives
---
title: "Assignment 1 - Binomial Model"
author: "Valery Pugacheva"
date: "January 30, 2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message=FALSE)
```


This is the solution to exercises 4 of Assignment 1. Let's first define the parameters:
```{r}
S0 <- 40
K <- 40
r <- .04
sigma <- .3
time_T <- 1/2
```

# Exercise 1: European Options

Let's first compute the limiting price of an European call/put (this is the celebrated Black-Scholes formula which you'll derive later in the course):
```{r}
d_1 <- (log(S0/K)+(r+sigma^2/2)*time_T) / (sigma*sqrt(time_T))
d_2 <- d_1 - sigma * sqrt(time_T)
bs_call_price <- S0*pnorm(d_1) - K*exp(-r*time_T)*pnorm(d_2)
bs_put_price <- K*exp(-r*time_T)*pnorm(-d_2) - S0*pnorm(-d_1)
bs_call_price
bs_put_price
```

To limit code duplication, it's convenient to define a function computing the no-arbitrage price:
```{r, results="hide"}
library(expm) # Raise matrix to a power

compute_no_arb_price <- function(payoff, N){
  delta_t <- time_T / N
  u <- exp(sigma*sqrt(delta_t))
  d <- exp(-sigma*sqrt(delta_t))
  q <- (exp(r*delta_t)-d) / (u-d)
  V <- rep(0, N+1)
  
  for (i in 1:(N+1)){
    V[i] <- payoff(S0 * u^(N+1-i) * d^(i-1))
    }
  Phi <- diag(q, N+1)
  Phi[col(Phi)==row(Phi)+1] <- 1-q
  Phi <- exp(-r*delta_t) * Phi
  return((Phi %^% N %*% V)[1])
}
```

Now can run it for puts and calls and plot the results. The horizontal line represents the limiting price computed in the beginning.

```{r, fig.height=3, echo=FALSE}
call_payoff <- function(S) pmax(S-K, 0)
put_payoff <- function(S) pmax(K-S, 0)

number_periods <- 50
call_prices <- matrix(nrow=2, ncol=number_periods)
call_prices[1,] <- 1:number_periods
put_prices <- matrix(nrow=2, ncol=number_periods)
put_prices[1,] <- 1:number_periods
for (N in 1:number_periods){
  call_prices[2,N] <- compute_no_arb_price(call_payoff, N)
  put_prices[2,N] <- compute_no_arb_price(put_payoff, N)
}

plot(call_prices[1,], call_prices[2,], type="l", xlab="# periods", 
     ylab="European call price")
abline(h=bs_call_price)

plot(put_prices[1,], put_prices[2,], type="l", xlab="# periods", 
     ylab="European put price")
```

If we use the proposed speed-up, we get the following:

```{r, fig.height=3, echo=FALSE}
number_periods <- 50
call_prices_speed <- matrix(nrow=2, ncol=number_periods)
call_prices_speed[1,] <- 1:number_periods
put_prices_speed <- matrix(nrow=2, ncol=number_periods)
put_prices_speed[1,] <- 1:number_periods
for (N in 1:number_periods){
  call_prices_speed[2,N] <- (compute_no_arb_price(call_payoff, N+1)+compute_no_arb_price(call_payoff, N))/2
  put_prices_speed[2,N] <- (compute_no_arb_price(put_payoff, N+1)+compute_no_arb_price(put_payoff, N))/2
}

plot(call_prices_speed[1,], call_prices_speed[2,], type="l", xlab="# periods", 
     ylab="European call price")
abline(h=bs_call_price)

plot(put_prices_speed[1,], put_prices_speed[2,], type="l", xlab="# periods", 
     ylab="European put price")
```

# Exercise 2: American Options

We need to modify our previous code:
```{r, results="hide"}
library(expm) # Raise matrix to a power

compute_american_no_arb_price <- function(payoff, N){
  delta_t <- time_T / N
  u <- exp(sigma*sqrt(delta_t))
  d <- exp(-sigma*sqrt(delta_t))
  q <- (exp(r*delta_t)-d) / (u-d)
  
  S <- matrix(0, nrow=N+1, ncol=N+1)
  for (i in 1:(N+1)){
    for (j in i:(N+1)){
      S[i,j] <- S0*u^(j-i)*d^(i-1)
    }
  }
  
  V <- payoff(S[,N+1])
  Phi <- diag(q, N+1)
  Phi[col(Phi)==row(Phi)+1] <- 1-q
  Phi <- exp(-r*delta_t) * Phi
  
  for (i in seq(N+1,2,-1)){
    V <- pmax(payoff(S[,i-1]), Phi %*% V)
    }
 
  return(V[1])
}
```

We again plot the results.

```{r, fig.height=3, echo=FALSE}
number_periods <- 10000
american_call_prices <- matrix(nrow=2, ncol=number_periods)
american_call_prices[1,] <- 1:number_periods
american_put_prices <- matrix(nrow=2, ncol=number_periods)
american_put_prices[1,] <- 1:number_periods
for (N in 1:number_periods){
  american_call_prices[2,N] <- compute_american_no_arb_price(call_payoff, N)
  american_put_prices[2,N] <- compute_american_no_arb_price(put_payoff, N)
}

plot(american_call_prices[1,], american_call_prices[2,], type="l", xlab="# periods", 
     ylab="American call price")
abline(h=bs_call_price)

plot(american_put_prices[1,], american_put_prices[2,], type="l", xlab="# periods", 
     ylab="American put price")
american_put_prices
```

Note that the prices of the American and the European call coincide. Indeed:
```{r}
all(round(american_call_prices[2,], 5) == round(call_prices[2,], 5))
```

If we use the proposed speed-up, we get the following:

```{r, fig.height=3, echo=FALSE}
number_periods <- 300
american_call_prices_speed <- matrix(nrow=2, ncol=number_periods)
american_call_prices_speed[1,] <- 1:number_periods
american_put_prices_speed <- matrix(nrow=2, ncol=number_periods)
american_put_prices_speed[1,] <- 1:number_periods
for (N in 1:number_periods){
  american_call_prices_speed[2,N] <- (compute_american_no_arb_price(call_payoff, N+1)+compute_american_no_arb_price(call_payoff, N))/2
  american_put_prices_speed[2,N] <- (compute_american_no_arb_price(put_payoff, N+1)+compute_american_no_arb_price(put_payoff, N))/2
}

plot(american_call_prices_speed[1,], american_call_prices_speed[2,], type="l", xlab="# periods", 
     ylab="American call price")
abline(h=bs_call_price)

plot(american_put_prices_speed[1,], american_put_prices_speed[2,], type="l", xlab="# periods", 
     ylab="American put price")
```
