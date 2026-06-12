B(1, 1) = V_g + l_e
B(2, 1) = 0
B(1, 2) = lambda
B(2, 2) = 0
B(1, 3) = 0
B(2, 3) = -U - alpha_grv*l_e
B(1, 4) = 0
B(2, 4) = -1
if (PRESENT(f)) then
f(1) = -lambda*y_T_1/(alpha_omg*omega_c**2)
f(2) = 0
endif
