function omega_approx = get_omega_approx(base, str, bc_index, t_indices, ks, varargin)

  % Use varargin to pass in flag for pulse charging

  if length(varargin) == 1
    pulse_flag = varargin{1};
  else
    pulse_flag = 0;
  end

  z_p = base.z_p;
  z = base.z;
  D_p = base.D_p;
  D = base.D;
  rho_s = base.rho_s;
  n = base.n;
  alpha = base.alpha;
  gamma = base.gamma;
  beta_m = base.beta_m;
  beta_v = base.beta_v;
  Da = base.Da;
  t_s = base.t_s;
  J_applied = base.J_applied;
  xi_p = base.xi_p;

  x = base.x;
  N = length(x);

  if pulse_flag == 1
    delta_t_on = base.delta_t_on;
    gamma_dc = base.gamma_dc;
  end

  % Assert key parameters
  parameters = str2double(split(str, '_'));
  assert(parameters(1) == rho_s)
  assert(parameters(2) == Da)
  assert(parameters(3) == J_applied)
  assert(parameters(4) == N)
  if pulse_flag == 1
    assert(parameters(5) == round(delta_t_on/t_s, 3, 'significant'))
    assert(parameters(6) == gamma_dc)
  end

  t = base.t;
  c_1 = base.c(end,:);  % at x = 1 for all t
  cx_1 = base.cx(end,:);  % at x = 1 for all t
  phix_1 = base.phix(end,:);  % at x = 1 for all t
  j_0_1 = base.j_0_1;
  eta_1 = base.eta_1;
  ct_1 = base.ct_1;

  alpha_1 = D - D_p;
  alpha_2 = z_p*D_p - z*D;

  N_t = length(t_indices);
  omega_approx = cell(N_t, 1);

  parfor i = 1:N_t

    ks_i = ks{i};

    t_temp = t;
    c_1_temp = c_1;
    cx_1_temp = cx_1;
    phix_1_temp = phix_1;
    j_0_1_temp = j_0_1;
    eta_1_temp = eta_1;
    ct_1_temp = ct_1;

    t_index = t_indices(i);

    t_i = t_temp(t_index);

    if t_i == 0

      beta_1 = 1 + (rho_s+abs(rho_s))/2;

      c_1_i = beta_1;
      cx_1_i = 0;
      phix_1_i = 0;
      j_0_1_i = Da*n*(xi_p*(beta_1-rho_s))^(1-alpha);
      eta_1_i = 0;
      ct_1_i = 0;

    else

      c_1_i = c_1_temp(t_index);
      cx_1_i = cx_1_temp(t_index);
      phix_1_i = phix_1_temp(t_index);
      j_0_1_i = j_0_1_temp(t_index);
      eta_1_i = eta_1_temp(t_index);
      ct_1_i = ct_1_temp(t_index);

    end

    % Prevent c_1_i from becoming too small because it appears in denominators of some terms
    if c_1_i < eps
      c_1_i = eps;
    end

    alpha_3_N = -alpha*exp(-alpha*n*eta_1_i) - (1-alpha)*exp((1-alpha)*n*eta_1_i);
    alpha_5_N = alpha_2*c_1_i - z_p*D_p*rho_s;
    xi_1 = ct_1_i./(z*c_1_i*D*ks_i);
    xi_2 = -(z*phix_1_i+ks_i)./(z*c_1_i*ks_i);
    G_1_hat = alpha_3_N*n*(-phix_1_i-gamma*ks_i.^2/n) + exp(-alpha*n*eta_1_i)*cx_1_i/(c_1_i-rho_s);
    G_2_hat = exp(-alpha*n*eta_1_i)/(c_1_i-rho_s);
    G_3_hat = -alpha_3_N*n;
    temp_1 = (alpha_1-alpha_5_N*xi_2).*ks_i - alpha_2*phix_1_i;
    temp_2 = beta_v*j_0_1_i*(G_1_hat-xi_1*G_3_hat) - beta_m*alpha_5_N*xi_1.*ks_i;
    temp_3 = beta_v*j_0_1_i*(G_2_hat+xi_2*G_3_hat);
    numerator = temp_1.*temp_2;
    denominator = temp_3 - beta_m*temp_1;
    omega_approx_i = beta_m*(numerator./denominator - alpha_5_N*xi_1.*ks_i);

    % Set omega = 0 when k = 0
    k_0 = ismembertol(ks_i, 0);
    omega_approx_i(k_0) = 0;

    omega_approx{i} = omega_approx_i;

  end

end
