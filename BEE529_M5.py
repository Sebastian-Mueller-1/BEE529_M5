# Total population, N.
N = 59864
# Initial number of infected and recovered individuals, I0 and R0.
I0 = 1; R0 = 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta = 0.33
gamma = 1./10 
zeta = 1/100
m = 7.7/1000/365 # From wikipedia
b = 18/1000/365 # From wikipedia
f = 0.011 # Guess about this disease?

print('r_0 is', beta/gamma)
# A grid of time points (in days)
simulation_days = 41
t = np.linspace(0, simulation_days, simulation_days)
print(t)

# The SIR model differential equations.
def SIRmodel(y, t):
    S = y[0]
    I = y[1]
    R = y[2]
    dSdt = -beta*S*I/(S+I+R) + R*zeta - m*S + b*(S+I+R)
    dIdt = beta*S*I/(S+I+R) - gamma*I - m*I - f*I
    dRdt = gamma*I - R*zeta - m*R
    return [dSdt, dIdt, dRdt]

# Initial conditions vector
y0 = [S0, I0, R0]
# Integrate the SIR equations over the time grid, t.
Y3 = odeint(SIRmodel, y0, t)
S3 = Y3[:,0]
I3 = Y3[:,1]
R3 = Y3[:,2]

print(S3)
# Initial number of infected and recovered individuals, I0 and R0.
I0 = I3[40]; R0 = R3[40]
# Everyone else, S0, is susceptible to infection initially.
S0 = S3[40]

# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta = 0.22
gamma = 1./10 
zeta = 1/100
m = 7.7/1000/365 # From wikipedia
b = 18/1000/365 # From wikipedia
f = .01 # Guess about this disease?
v = .01 # vaccination coefficient assuming 1% of any given amount of individuals get vaccinated 

# A grid of time points (in days)
simulation_days = 400
t2 = np.linspace(42, simulation_days, simulation_days)

# The SIR model differential equations.
def SIRmodel(y, t):
    S = y[0]
    I = y[1]
    R = y[2]
    
    dSdt = -beta*S*I/(S+I+R) + R*zeta - m*S + b*(S+I+R) - v*S
    dIdt = beta*S*I/(S+I+R) - gamma*I - m*I - f*I
    dRdt = gamma*I - R*zeta - m*R + v*S
    
    return [dSdt, dIdt, dRdt]

# Initial conditions vector
y0 = [S0, I0, R0]
# Integrate the SIR equations over the time grid, t.
Y4 = odeint(SIRmodel, y0, t2)
S4 = Y4[:,0]
I4 = Y4[:,1]
R4 = Y4[:,2]

t_final = np.concatenate([t,t2])
s_final = np.concatenate([S3,S4])
i_final = np.concatenate([I3,I4])
r_final = np.concatenate([R3,R4])

plt.figure(1)
plt.plot(t_final, s_final/1000, 'b', label='Susceptible')
plt.plot(t_final, i_final/1000, 'r', label='Infected')
plt.plot(t_final, r_final/1000, 'g', label='Recovered with immunity')
plt.axvline(x=42, color = "black", linestyle= 'dashed', label = 'state change', alpha= .3)
plt.xlabel('Time /days')
plt.ylabel('Number (1000s)')
plt.legend(loc='upper right')
plt.show()