
if (trans_) then
B(1, 1) = 1
B(2, 1) = -1
B(3, 1) = 0
B(4, 1) = 0
B(1, 2) = U*alpha_gbc
B(2, 2) = 0
B(3, 2) = -U - l_e
B(4, 2) = -1
else
B(1, 1) = 1
B(2, 1) = U*alpha_gbc
B(1, 2) = -1
B(2, 2) = 0
B(1, 3) = 0
B(2, 3) = -U - l_e
B(1, 4) = 0
B(2, 4) = -1
endif
if (PRESENT(s)) then
s(1) = 0
s(2) = 0
endif
