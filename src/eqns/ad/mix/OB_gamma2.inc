
if (trans_) then
B(1, 1) = V_g + l_e
B(2, 1) = lambda/(alpha_omg*c_1*omega_c**2)
B(3, 1) = -lambda/(alpha_omg*c_1*omega_c**2)
B(4, 1) = 0
B(1, 2) = 0
B(2, 2) = 0
B(3, 2) = -U - l_e
B(4, 2) = -1
else
B(1, 1) = V_g + l_e
B(2, 1) = 0
B(1, 2) = lambda/(alpha_omg*c_1*omega_c**2)
B(2, 2) = 0
B(1, 3) = -lambda/(alpha_omg*c_1*omega_c**2)
B(2, 3) = -U - l_e
B(1, 4) = 0
B(2, 4) = -1
endif
if (PRESENT(s)) then
s(1) = -lambda*y_T_1/(alpha_omg*omega_c**2)
s(2) = 0
endif
