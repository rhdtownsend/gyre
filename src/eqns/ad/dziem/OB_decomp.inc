
if (trans_) then
B(1, 1) = a_11 - chi
B(2, 1) = a_12
B(3, 1) = -G_1*alpha_grv - a_12*alpha_grv
B(4, 1) = G_2*alpha_grv
B(1, 2) = U*alpha_gbc
B(2, 2) = 0
B(3, 2) = l_e + 1
B(4, 2) = 1
else
B(1, 1) = a_11 - chi
B(2, 1) = U*alpha_gbc
B(1, 2) = a_12
B(2, 2) = 0
B(1, 3) = -G_1*alpha_grv - a_12*alpha_grv
B(2, 3) = l_e + 1
B(1, 4) = G_2*alpha_grv
B(2, 4) = 1
endif
if (PRESENT(s)) then
s(1) = c_1*y_T_1*(G_1 - G_2*l_e)
s(2) = 0
endif
