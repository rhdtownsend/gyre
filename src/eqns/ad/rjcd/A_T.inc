A_T(1, 1) = invx*(V_g - l_i - 1)
A_T(2, 1) = invx*om2*(-V_g*alpha_gam + eta)
A_T(3, 1) = invx*(-V_g*alpha_gam - alpha_grv*eta + eta)
A_T(4, 1) = 0
A_T(1, 2) = invx*(-As_iso + om2)/om2
A_T(2, 2) = invx*(-om2*(3 - U) + om2*(As - U - l_i + 3))/om2
A_T(3, 2) = invx*((alpha_grv*(1 - U) - alpha_grv*(-U - l_i + 3))/om2 + ( &
      As - U - alpha_grv*(1 - U) - l_i + 3)/om2)
A_T(4, 2) = 0
A_T(1, 3) = 0
A_T(2, 3) = 0
A_T(3, 3) = invx*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A_T(4, 3) = alpha_grv*invx
A_T(1, 4) = -As*U*alpha_grv*invx
A_T(2, 4) = -U*V_g*alpha_grv*invx*om2
A_T(3, 4) = invx*(-U*V_g*alpha_grv - U*(-As - U - V_g + 3) + alpha_grv* &
      lambda - alpha_grv*(1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 &
      - U) - alpha_grv*(-U - l_i + 3)))
A_T(4, 4) = invx*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
