
if (trans_) then
A(1, 1) = del_x*(-U - l_e + 3)
A(2, 1) = del_x
A(1, 2) = del_x*(-As*U - U*V_g + lambda)
A(2, 2) = del_x*(-U - l_e + 2)
else
A(1, 1) = del_x*(-U - l_e + 3)
A(2, 1) = del_x*(-As*U - U*V_g + lambda)
A(1, 2) = del_x
A(2, 2) = del_x*(-U - l_e + 2)
endif
if (PRESENT(r)) then
r(1) = 0
r(2) = c_1*del_x*y_T_1*(-As*U - U*V_g)
endif
