
if (trans_) then
B(1, 1) = a_11 + a_12 - chi
B(2, 1) = a_12/V_2
else
B(1, 1) = a_11 + a_12 - chi
B(1, 2) = a_12/V_2
endif
if (PRESENT(s)) then
s(1) = 0
endif
