A_t(1, 1) = -del_x
A_t(2, 1) = -V_g*del_x/V_2
A_t(1, 2) = del_x*(V_2*(om2 + 3) + V_2)
A_t(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U + 3)/V_2) + V_g)
if (PRESENT(f)) then
f(1) = 0
f(2) = 0
endif
