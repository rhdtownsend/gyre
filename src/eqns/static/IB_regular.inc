
if (trans_) then
B(1, 1) = l_e
B(2, 1) = -1
else
B(1, 1) = l_e
B(1, 2) = -1
endif
if (PRESENT(s)) then
s(1) = 0
endif
