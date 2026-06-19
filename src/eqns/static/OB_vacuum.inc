
if (trans_) then
B(1, 1) = -U + l_e + 1
B(2, 1) = 1
else
B(1, 1) = -U + l_e + 1
B(1, 2) = 1
endif
if (PRESENT(s)) then
s(1) = -U*y_T_1
endif
