
if (trans_) then
B(1, 1) = a_11 + a_12 - chi
B(2, 1) = a_12/V_2
B(3, 1) = -G_1*alpha_grv
B(4, 1) = G_2*alpha_grv
B(5, 1) = 0
B(6, 1) = 0
B(1, 2) = U*alpha_gbc
B(2, 2) = 0
B(3, 2) = l_e + 1
B(4, 2) = 1
B(5, 2) = 0
B(6, 2) = 0
B(1, 3) = 2
B(2, 3) = 4*V*nabla_ad/V_2
B(3, 3) = 0
B(4, 3) = 0
B(5, 3) = 4*f_rht
B(6, 3) = -1
else
B(1, 1) = a_11 + a_12 - chi
B(2, 1) = U*alpha_gbc
B(3, 1) = 2
B(1, 2) = a_12/V_2
B(2, 2) = 0
B(3, 2) = 4*V*nabla_ad/V_2
B(1, 3) = -G_1*alpha_grv
B(2, 3) = l_e + 1
B(3, 3) = 0
B(1, 4) = G_2*alpha_grv
B(2, 4) = 1
B(3, 4) = 0
B(1, 5) = 0
B(2, 5) = 0
B(3, 5) = 4*f_rht
B(1, 6) = 0
B(2, 6) = 0
B(3, 6) = -1
endif
if (PRESENT(s)) then
s(1) = c_1*y_T_1*(G_1 - G_2*l_e)
s(2) = 0
s(3) = 0
endif
