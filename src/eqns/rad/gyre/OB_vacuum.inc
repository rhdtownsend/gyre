
if (trans_) then
B(1, 1) = 1
B(2, 1) = -1
else
B(1, 1) = 1
B(1, 2) = -1
endif
if (PRESENT(s)) then
s(1) = 0
endif
