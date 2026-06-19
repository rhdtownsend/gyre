
if (trans_) then
C(1, 1) = 1
C(2, 1) = 0
C(1, 2) = 0
C(2, 2) = U/V_2
else
C(1, 1) = 1
C(2, 1) = 0
C(1, 2) = 0
C(2, 2) = U/V_2
endif
if (PRESENT(t)) then
t(1) = 0
t(2) = 0
endif
