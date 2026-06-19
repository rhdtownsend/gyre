
if (trans_) then
B(1, 1) = -b_11 - b_12 + chi
B(2, 1) = -b_12/V_2
else
B(1, 1) = -b_11 - b_12 + chi
B(1, 2) = -b_12/V_2
endif
if (PRESENT(s)) then
s(1) = 0
endif
