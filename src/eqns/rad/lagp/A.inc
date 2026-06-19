
if (trans_) then
A(1, 1) = -del_x
A(2, 1) = -V_g*del_x/V_2
A(1, 2) = del_x*(V_2*(om2 + 3) + V_2)
A(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U + 3)/V_2) + V_g)
else
A(1, 1) = -del_x
A(2, 1) = del_x*(V_2*(om2 + 3) + V_2)
A(1, 2) = -V_g*del_x/V_2
A(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U + 3)/V_2) + V_g)
endif
if (PRESENT(r)) then
r(1) = 0
r(2) = 0
endif
