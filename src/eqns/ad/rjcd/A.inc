A(1, 1) = invx*(V_g - l_i - 1)
A(2, 1) = invx*(-As_iso + om2)/om2
A(3, 1) = 0
A(4, 1) = -As*U*alpha_grv*invx
A(1, 2) = invx*om2*(-V_g*alpha_gam + eta)
A(2, 2) = invx*(-om2*(3 - U) + om2*(As - U - l_i + 3))/om2
A(3, 2) = 0
A(4, 2) = -U*V_g*alpha_grv*invx*om2
A(1, 3) = invx*(-V_g*alpha_gam - alpha_grv*eta + eta)
A(2, 3) = invx*((alpha_grv*(1 - U) - alpha_grv*(-U - l_i + 3))/om2 + (As &
      - U - alpha_grv*(1 - U) - l_i + 3)/om2)
A(3, 3) = invx*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A(4, 3) = invx*(-U*V_g*alpha_grv - U*(-As - U - V_g + 3) + alpha_grv* &
      lambda - alpha_grv*(1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 &
      - U) - alpha_grv*(-U - l_i + 3)))
A(1, 4) = 0
A(2, 4) = 0
A(3, 4) = alpha_grv*invx
A(4, 4) = invx*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
