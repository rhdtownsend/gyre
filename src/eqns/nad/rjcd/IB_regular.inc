
if (trans_) then
B(1, 1) = alpha_omg*c_1*omega_c**2
B(2, 1) = -alpha_omg*c_1*l_i*omega_c**2
B(3, 1) = alpha_grv*l_i - l_i
B(4, 1) = 0
B(5, 1) = 0
B(6, 1) = 0
B(1, 2) = 0
B(2, 2) = 0
B(3, 2) = -alpha_grv*l_i + 4*alpha_grv - 2
B(4, 2) = 2*alpha_grv - 1
B(5, 2) = 0
B(6, 2) = 0
B(1, 3) = 0
B(2, 3) = 0
B(3, 3) = 0
B(4, 3) = 0
B(5, 3) = 1
B(6, 3) = 0
else
B(1, 1) = alpha_omg*c_1*omega_c**2
B(2, 1) = 0
B(3, 1) = 0
B(1, 2) = -alpha_omg*c_1*l_i*omega_c**2
B(2, 2) = 0
B(3, 2) = 0
B(1, 3) = alpha_grv*l_i - l_i
B(2, 3) = -alpha_grv*l_i + 4*alpha_grv - 2
B(3, 3) = 0
B(1, 4) = 0
B(2, 4) = 2*alpha_grv - 1
B(3, 4) = 0
B(1, 5) = 0
B(2, 5) = 0
B(3, 5) = 1
B(1, 6) = 0
B(2, 6) = 0
B(3, 6) = 0
endif
if (PRESENT(s)) then
s(1) = c_1*l_i*y_T_1
s(2) = 0
s(3) = 0
endif
