
if (trans_) then
C(1, 1) = 1
C(2, 1) = 0
C(1, 2) = -U
C(2, 2) = U
else
C(1, 1) = 1
C(2, 1) = -U
C(1, 2) = 0
C(2, 2) = U
endif
if (PRESENT(F)) then
F(1) = 0
F(2) = 0
endif
