
if (trans_) then
B(1, 1) = 0
B(2, 1) = -1/V_2
else
B(1, 1) = 0
B(1, 2) = -1/V_2
endif
if (PRESENT(s)) then
s(1) = 0
endif
