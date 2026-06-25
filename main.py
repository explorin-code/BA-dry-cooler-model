import matplotlib.pyplot as plt
from solvers import solve_it_LMTD

# Run the solver
k, Q, T_c_out, T_a_out, hist_hot, hist_cold = solve_it_LMTD()

print(f"Solver converged! k = {k:.2f} W/m²K, Q = {Q/1000:.2f} kW")

# Plotting the oscillation
plt.figure(figsize=(10, 5))
plt.plot(hist_hot, label='Error dT_hot [K]', marker='o')
plt.plot(hist_cold, label='Error dT_cold [K]', marker='o')
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.title("Solver Convergence Oscillation")
plt.xlabel("Iteration")
plt.ylabel("Difference (Current - Previous) [K]")
plt.legend()
plt.grid(True)
plt.show()