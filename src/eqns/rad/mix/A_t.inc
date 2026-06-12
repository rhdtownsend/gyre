A_t(1, 1) = del_x*(V_g - 1)
A_t(2, 1) = -V_g*del_x
A_t(1, 2) = del_x*(-As + U + om2)
A_t(2, 2) = del_x*(As - U + 3)
if (PRESENT(f)) then
f(1) = 0
f(2) = 0
endif
