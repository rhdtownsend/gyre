
if (trans_) then
A(1, 1) = del_x*(V_g - 1)
A(2, 1) = -V_g*del_x
A(1, 2) = del_x*(-As + U + om2)
A(2, 2) = del_x*(As - U + 3)
else
A(1, 1) = del_x*(V_g - 1)
A(2, 1) = del_x*(-As + U + om2)
A(1, 2) = -V_g*del_x
A(2, 2) = del_x*(As - U + 3)
endif
if (PRESENT(r)) then
r(1) = 0
r(2) = 0
endif
