A_T(1, 1) = invx*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A_T(2, 1) = invx*(-V_g*alpha_gam + eta)/V_2
A_T(3, 1) = alpha_grv*eta*invx
A_T(4, 1) = 0
A_T(1, 2) = invx*(-V_2*(-V_g*alpha_gam + V_g + eta - l_i - 1) + V_2*(As &
      - As_iso - U - l_i + om2 + 3))
A_T(2, 2) = invx*(V_2*(dV_2/V_2 + (As - U - l_i + 3)/V_2) + V_g* &
      alpha_gam - eta)
A_T(3, 2) = -V_2*alpha_grv*eta*invx
A_T(4, 2) = -V_2*alpha_grv*invx
A_T(1, 3) = 0
A_T(2, 3) = 0
A_T(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A_T(4, 3) = alpha_grv*invx
A_T(1, 4) = invx*(As*U*alpha_grv + U*V_g*alpha_grv)
A_T(2, 4) = U*V_g*alpha_grv*invx/V_2
A_T(3, 4) = alpha_grv*invx*lambda
A_T(4, 4) = alpha_grv*invx*(-U - l_i + 2)
