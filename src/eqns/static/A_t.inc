A_t(1, 1) = del_x*(-U - l_e + 3)
A_t(2, 1) = del_x
A_t(1, 2) = del_x*(-As*U - U*V_g + lambda)
A_t(2, 2) = del_x*(-U - l_e + 2)
if (PRESENT(f)) then
f(1) = 0
f(2) = c_1*del_x*y_T_1*(-As*U - U*V_g)
endif
