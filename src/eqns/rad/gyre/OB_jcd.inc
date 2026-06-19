
if (trans_) then
B(1, 1) = -b_11 + chi
B(2, 1) = -b_12
else
B(1, 1) = -b_11 + chi
B(1, 2) = -b_12
endif
if (PRESENT(s)) then
s(1) = 0
endif
