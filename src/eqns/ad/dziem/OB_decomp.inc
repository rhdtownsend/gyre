B(1, 1) = a_11 - chi
B(2, 1) = U*alpha_gbc*alpha_grv
B(1, 2) = a_12
B(2, 2) = 0
B(1, 3) = -G_1*alpha_grv - a_12
B(2, 3) = alpha_grv*l_e + 1
B(1, 4) = G_2*alpha_grv
B(2, 4) = alpha_grv
if (PRESENT(f)) then
f(1) = c_1*y_T_1*(G_1 - G_2*l_e)
f(2) = 0
endif
