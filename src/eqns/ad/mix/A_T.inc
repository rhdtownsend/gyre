A_T(1, 1) = invx*(V_g - l_i - 1)
A_T(2, 1) = invx*(-V_g*alpha_gam + eta)
A_T(3, 1) = -alpha_grv*eta*invx
A_T(4, 1) = 0
A_T(1, 2) = invx*(-As_iso + om2)
A_T(2, 2) = invx*(As - U - l_i + 3)
A_T(3, 2) = -alpha_grv*invx*(1 - U)
A_T(4, 2) = alpha_grv*invx
A_T(1, 3) = 0
A_T(2, 3) = 0
A_T(3, 3) = invx*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A_T(4, 3) = alpha_grv*invx
A_T(1, 4) = -As*U*alpha_grv*invx
A_T(2, 4) = -U*V_g*alpha_grv*invx
A_T(3, 4) = invx*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv* &
      (1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv*( &
      -U - l_i + 3)))
A_T(4, 4) = invx*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
