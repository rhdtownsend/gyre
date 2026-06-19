
if (trans_) then
B(1, 1) = 1
B(2, 1) = 0
B(3, 1) = 0
B(4, 1) = 0
B(1, 2) = U*alpha_gbc*alpha_grv
B(2, 2) = 0
B(3, 2) = alpha_grv*l_e + 1
B(4, 2) = alpha_grv
else
B(1, 1) = 1
B(2, 1) = U*alpha_gbc*alpha_grv
B(1, 2) = 0
B(2, 2) = 0
B(1, 3) = 0
B(2, 3) = alpha_grv*l_e + 1
B(1, 4) = 0
B(2, 4) = alpha_grv
endif
if (PRESENT(s)) then
s(1) = 0
s(2) = 0
endif
