# Set parameters
S0 <- 110
T <- 1
r <- 0.01
sigma <- 0.2
K <- 100
m <- 12 # Number of steps in one simulation
deltaT <- T/m # length of one step


# Simulate stock price
Stock <- function(st,r,sigma,deltaT){
  st*exp((r-1/2*sigma^2)*deltaT+sigma*sqrt(deltaT)*rnorm(1))
}

# Wrap the codes in a function MC_AsianCall
MC_AsianCall <- function(N){
########## STOCK PRICE PATHS ####################
# Simulate paths of stock price
Stock_Path <- function(x){
  stock_path <- numeric(m+1)
  stock_path[1] <- S0
  for (y in 1:m){
    stock_path[1+y] <- Stock(st=stock_path[y],r,sigma,deltaT)
  }
  return(stock_path)
}

# Use matrix Stock_Path_Matrix to store all of N simulations result
Stock_Path_Matrix <- sapply(1:N,Stock_Path)
colnames(Stock_Path_Matrix) <- paste0("Path", 1:ncol(Stock_Path_Matrix))
rownames(Stock_Path_Matrix) <- paste0("S", 0:(nrow(Stock_Path_Matrix)-1))

########## ARITHMETIC ASIAN OPTION ####################
# Caculate the price of each simulation(path)
Arith_AsianOption_Price <- colMeans(Stock_Path_Matrix[-1,])-K  # Get rid of S0
Arith_AsianOption_Price[Arith_AsianOption_Price<0] <- 0  # Set negative payoff to 0
Arith_AsianOption_Price <- exp(-r*T)*Arith_AsianOption_Price # Discount

########## GEOMETRIC ASIAN OPTION ####################
# Caculate theoretical price of geometric asian option based on Black-Scholes formula
# Notice t(i)=i*deltaT=i*T/m
t <- numeric(m)
for (i in 1:m){t[i]<-i*deltaT}
T_bar <- mean(t)
for (i in 1:m){t[i]<-(2*i-1)*(m+1-i)*deltaT}
sigma_bar <- sqrt(sigma^2/(m^2*T_bar)*sum(t))
delta <- 0.5*sigma^2-0.5*sigma_bar^2
d1 <- (log(S0/K)+(r-delta+0.5*sigma_bar^2)*T_bar)/(sigma_bar*sqrt(T_bar))
BS_Geo_AsianOption_Price <- exp(-delta*T_bar-r*(T-T_bar))*S0*pnorm(d1)-exp(-r*T)*K*pnorm(d1-sigma_bar*sqrt(T_bar))

# Caculate the price of each simulation(path)
Geo_AsianOption_Price <- (apply(Stock_Path_Matrix[-1,],2,prod))^(1/m)-K
Geo_AsianOption_Price[Geo_AsianOption_Price<0] <- 0  # Set negative payoff to 0
Geo_AsianOption_Price <- exp(-r*T)*Geo_AsianOption_Price # Discount

########## CACULATE RESULT ####################
b_star <- cov(x=Geo_AsianOption_Price,y=Arith_AsianOption_Price)/var(x=Geo_AsianOption_Price)
Var_Reduced_Arith_AsianOption_Price <- Arith_AsianOption_Price-b_star*(Geo_AsianOption_Price-BS_Geo_AsianOption_Price)

error1 <- sd(Arith_AsianOption_Price)/sqrt(N)
error2 <- sd(Var_Reduced_Arith_AsianOption_Price)/sqrt(N)
Correlation <- cor(Geo_AsianOption_Price,Arith_AsianOption_Price)
c("Number of replications"=N,"1st Algo Price (Arithmetic)"=mean(Arith_AsianOption_Price)
  ,"2nd Algo Price (Variance Reduction)"=mean(Var_Reduced_Arith_AsianOption_Price)
  ,"Standard Error of 1st Algo"=error1,"Standard Error of 2nd Algo"=error2
  ,"Correlation"=Correlation)
}

########## DISPLAY RESULT ####################
result <- lapply(10^(2:6),MC_AsianCall)
result <- do.call(rbind,result)
result
write.csv(result,"C:/raymond/NYU Tandon/Stochastic Calculus and Option pricing FRE6233/HW/result.csv",row.names = F) 
# Conclusion: With the increment of number of simulations, the option prices from both
# 1st and 2nd algorithms tend to be stable, and the standard errors of both algos 
# become smaller. Notice that the standard error of 2nd algo (variance reduced)
# is always smaller than the one of 1st algo (arithmetic Asian call). And the
# correlation between arithmetic and geometric Asian call is always higher than 0.999

# Note: My simulation doesn't include S0=110, so there are only 12 steps instead of 13
# steps in one simulation. The option price will change a little if there are 13 steps.