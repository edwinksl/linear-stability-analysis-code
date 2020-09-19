function copper(save_flag)

  close all

  voltammetry_path = '/home/edwinksl/git/voltammetry/';
  ss_dir = 'code/1d-steady-state/';
  copper_dir = 'datasets/copper/';

  % Add path
  addpath([voltammetry_path, ss_dir]);

  % Get parameters CN_2(-) and CN_2(+)
  params_CN_b_m = system_parameters(4, 'm');
  params_CN_b_p = system_parameters(4, 'p');

  % Get relevant parameters
  z_p = params_CN_b_m.z_p;
  z_m = params_CN_b_m.z_m;
  beta_m = params_CN_b_m.beta_m;
  Omega = params_CN_b_m.Omega_m;  % [=] m^3
  L_x = params_CN_b_m.L;  % [=] m
  r_p = params_CN_b_m.r_p;  % [=] m

  D_p_opt_CN_b_m = params_CN_b_m.D_p_0_tilde_opt;
  D_m_opt_CN_b_m = params_CN_b_m.D_m_0_tilde_opt;
  rho_s_opt_CN_b_m = params_CN_b_m.rho_s_tilde_opt;
  J_0_ref_tilde_opt_CN_b_m = params_CN_b_m.J_0_ref_tilde_opt;
  Da_opt_CN_b_m = params_CN_b_m.Da_opt;
  alpha_1_opt_CN_b_m = params_CN_b_m.alpha_1_opt;
  I_lim_opt_CN_b_m = params_CN_b_m.I_lim_opt;
  dt_opt_CN_b_m = params_CN_b_m.diffusion_time_opt;
  omega_scale_opt_CN_b_m = params_CN_b_m.omega_scale_opt;

  rho_s_opt_CN_b_m = round(rho_s_opt_CN_b_m, 2);  % magnitude of original value is too small, so make it slightly bigger by rounding

  D_p_opt_CN_b_p = params_CN_b_p.D_p_0_tilde_opt;
  D_m_opt_CN_b_p = params_CN_b_p.D_m_0_tilde_opt;
  rho_s_opt_CN_b_p = params_CN_b_p.rho_s_tilde_opt;
  J_0_ref_tilde_opt_CN_b_p = params_CN_b_p.J_0_ref_tilde_opt;
  Da_opt_CN_b_p = params_CN_b_p.Da_opt;
  alpha_1_opt_CN_b_p = params_CN_b_p.alpha_1_opt;
  I_lim_opt_CN_b_p = params_CN_b_p.I_lim_opt;
  dt_opt_CN_b_p = params_CN_b_p.diffusion_time_opt;
  omega_scale_opt_CN_b_p = params_CN_b_p.omega_scale_opt;

  I_max_opt_CN_b_p = get_I_max(z_p, z_m, rho_s_opt_CN_b_p, 3);
  I_max_opt_CN_b_p_dim = I_max_opt_CN_b_p*I_lim_opt_CN_b_p;

  % Assign BV parameters to their optimal values
  params_CN_b_m.J_0_ref_tilde = J_0_ref_tilde_opt_CN_b_m;
  params_CN_b_m.alpha_1 = alpha_1_opt_CN_b_m;

  params_CN_b_p.J_0_ref_tilde = J_0_ref_tilde_opt_CN_b_p;
  params_CN_b_p.alpha_1 = alpha_1_opt_CN_b_p;

  % Print parameters
  fprintf('CN_b_m\n')
  fprintf('D_p = %#.3g\n', D_p_opt_CN_b_m)
  fprintf('D_m = %#.3g\n', D_m_opt_CN_b_m)
  fprintf('rho_s = %#.3g\n', rho_s_opt_CN_b_m)
  fprintf('J_0_ref_tilde = %#.3g\n', J_0_ref_tilde_opt_CN_b_m)
  fprintf('Da = %#.3g\n', Da_opt_CN_b_m)
  fprintf('alpha_1 = %#.3g\n', alpha_1_opt_CN_b_m)
  fprintf('I_lim = %#.3g mA\n', I_lim_opt_CN_b_m*1e3)
  fprintf('Diffusion time = %#.3g s\n', dt_opt_CN_b_m)
  fprintf('Omega scale = %#.3g s^-1\n', omega_scale_opt_CN_b_m)

  fprintf('\n')

  fprintf('CN_b_p\n')
  fprintf('D_p = %#.3g\n', D_p_opt_CN_b_p)
  fprintf('D_m = %#.3g\n', D_m_opt_CN_b_p)
  fprintf('rho_s = %#.3g\n', rho_s_opt_CN_b_p)
  fprintf('J_0_ref_tilde = %#.3g\n', J_0_ref_tilde_opt_CN_b_p)
  fprintf('Da = %#.3g\n', Da_opt_CN_b_p)
  fprintf('alpha_1 = %#.3g\n', alpha_1_opt_CN_b_p)
  fprintf('I_lim = %#.3g mA\n', I_lim_opt_CN_b_p*1e3)
  fprintf('I_max = %#.3g mA\n', I_max_opt_CN_b_p_dim*1e3)
  fprintf('Diffusion time = %#.3g s\n', dt_opt_CN_b_p)
  fprintf('Omega scale = %#.3g s^-1\n', omega_scale_opt_CN_b_p)

  fprintf('\n')

  % Physical constants
  k_B = 1.38064852e-23;  % [=] m^2 kg s^-2 K^-1

  % Dimensional LSA parameters
  T = 298;  % [=] K
  n = 2;
  gamma_dim = 1.85;  % [=] J m^-2
  d_p = 2*r_p;  % [=] m
  d_p_lb = 200e-9;  % [=] m
  d_p_ub = 300e-9;  % [=] m
  h_c = 2*d_p;  % [=] m
  h_c_lb = 2*d_p_lb;  % [=] m
  h_c_ub = 2*d_p_ub;  % [=] m
  L_y = 2*6.5e-3;  % [=] m; 2*6.5 mm
  L_z = 2*6.5e-3;  % [=] m; 2*6.5 mm
  I_dims = [15e-3, 20e-3, 25e-3];  % [=] A

  data.h_c = h_c;
  data.h_c_lb = h_c_lb;
  data.h_c_ub = h_c_ub;
  data.I_dims = I_dims;

  % Dimensionless LSA parameters
  gamma = gamma_dim/(L_x*k_B*T/Omega);

  fprintf('gamma = %#.3g\n', gamma)

  fprintf('\n')

  beta_D_opt_CN_b_m = get_beta_D(D_p_opt_CN_b_m, D_m_opt_CN_b_m, z_p, z_m);
  beta_D_opt_CN_b_p = get_beta_D(D_p_opt_CN_b_p, D_m_opt_CN_b_p, z_p, z_m);

  beta_v_opt_CN_b_m = beta_m/beta_D_opt_CN_b_m;
  beta_v_opt_CN_b_p = beta_m/beta_D_opt_CN_b_p;

  Is_CN_b_m = I_dims/I_lim_opt_CN_b_m;
  Is_CN_b_p = I_dims/I_lim_opt_CN_b_p;

  N_I = length(I_dims);
  k_c_dims = NaN(2, N_I);
  lambda_c_dims = NaN(2, N_I);
  k_max_dims = NaN(2, N_I);
  lambda_max_dims = NaN(2, N_I);
  omega_max_dims = NaN(2, N_I);

  for i = 1:2

    if i == 1
      params = params_CN_b_m;
      D_p_opt = D_p_opt_CN_b_m;
      D_m_opt = D_m_opt_CN_b_m;
      rho_s_opt = rho_s_opt_CN_b_m;
      alpha_1_opt = alpha_1_opt_CN_b_m;
      omega_scale_opt = omega_scale_opt_CN_b_m;
      beta_v_opt = beta_v_opt_CN_b_m;
      Is = Is_CN_b_m;
    elseif i == 2
      params = params_CN_b_p;
      D_p_opt = D_p_opt_CN_b_p;
      D_m_opt = D_m_opt_CN_b_p;
      rho_s_opt = rho_s_opt_CN_b_p;
      alpha_1_opt = alpha_1_opt_CN_b_p;
      omega_scale_opt = omega_scale_opt_CN_b_p;
      beta_v_opt = beta_v_opt_CN_b_p;
      Is = Is_CN_b_p;
    end

    for j = 1:N_I

      I = Is(j);

      fprintf('rho_s = %#.3g\n', rho_s_opt)
      fprintf('I = %#.3g\n', I)

      fprintf('\n')

      if i == 2 && I > I_max_opt_CN_b_p

        k_c_approx = Inf;
        k_max_approx = Inf;
        omega_max_approx = Inf;

        fprintf('I > I_max for rho_s > 0\n')
        fprintf('rho_s = %#.3g\n', rho_s_opt)
        fprintf('I = %#.3g\n', I)
        fprintf('I_max = %#.3g\n', I_max_opt_CN_b_p)
        fprintf('\n')

      else

        % Compute base state variables at cathode (x = 1)
        [c_0, c_0_x, ~, phi_0, phi_0_x, ~, V_0] = get_ss_derivatives(1, z_p, z_m, I, 1, rho_s_opt, 3, params);
        phi_e_c = -V_0;
        [~, ~, eta_0, ~, j_0_0] = bv_ss(c_0, phi_0, rho_s_opt, phi_e_c, params);

        % Compute k_c_approx
        k_c_params.rho_s = rho_s_opt;
        k_c_params.n = n;
        k_c_params.alpha = alpha_1_opt/2;
        k_c_params.gamma = 2*gamma;
        k_c_params.c_1 = c_0;
        k_c_params.cx_1 = c_0_x;
        k_c_params.phix_1 = phi_0_x;
        k_c_params.eta_1 = eta_0;

        k_c_approx = get_k_c_approx(k_c_params);

        % Compute max_approx
        k_max_params.z_p = z_p;
        k_max_params.z = z_m;
        k_max_params.D_p = D_p_opt;
        k_max_params.D = D_m_opt;
        k_max_params.rho_s = rho_s_opt;
        k_max_params.n = n;
        k_max_params.alpha = alpha_1_opt/2;
        k_max_params.gamma = 2*gamma;
        k_max_params.beta_m = beta_m;
        k_max_params.beta_v = beta_v_opt;

        k_max_params.c_1 = c_0;
        k_max_params.cx_1 = c_0_x;
        k_max_params.phix_1 = phi_0_x;
        k_max_params.j_0_1 = j_0_0;
        k_max_params.eta_1 = eta_0;
        k_max_params.ct_1 = 0;

        [k_max_approx, omega_max_approx] = get_max_approx(k_max_params, k_c_approx);

      end

      k_c_approx_dim = k_c_approx/L_x;
      lambda_c_approx_dim = 2*pi/k_c_approx_dim;

      fprintf('k_c_approx = %#.3g\n', k_c_approx)
      fprintf('k_c_approx_dim = %#.3g m^-1\n', k_c_approx_dim)
      fprintf('h_c = %#.3g um\n', h_c*1e6)
      fprintf('h_c_lb = %#.3g um\n', h_c_lb*1e6)
      fprintf('h_c_ub = %#.3g um\n', h_c_ub*1e6)
      fprintf('lambda_c_approx_dim = %#.3g um\n', lambda_c_approx_dim*1e6)

      fprintf('\n')

      k_max_approx_dim = k_max_approx/L_x;
      lambda_max_approx_dim = 2*pi/k_max_approx_dim;
      omega_max_approx_dim = omega_max_approx*omega_scale_opt;

      fprintf('k_max_approx = %#.3g\n', k_max_approx)
      fprintf('k_max_approx_dim = %#.3g m^-1\n', k_max_approx_dim)
      fprintf('lambda_max_approx_dim = %#.3g um\n', lambda_max_approx_dim*1e6)
      fprintf('omega_max_approx = %#.3g\n', omega_max_approx)
      fprintf('omega_max_approx_dim = %#.2e s^-1\n', omega_max_approx_dim)

      fprintf('\n')

      k_c_dims(i,j) = k_c_approx_dim;
      lambda_c_dims(i,j) = lambda_c_approx_dim;
      k_max_dims(i,j) = k_max_approx_dim;
      lambda_max_dims(i,j) = lambda_max_approx_dim;
      omega_max_dims(i,j) = omega_max_approx_dim;

    end

  end

  data.k_c_dims = k_c_dims;
  data.lambda_c_dims = lambda_c_dims;
  data.k_max_dims = k_max_dims;
  data.lambda_max_dims = lambda_max_dims;
  data.omega_max_dims = omega_max_dims;

  % Save data
  if save_flag == 1
    data_filename = 'copper.mat';
    save([copper_dir, data_filename], 'data')
  end

end

function k_c_approx = get_k_c_approx(params)

  rho_s = params.rho_s;
  n = params.n;
  alpha = params.alpha;
  gamma = params.gamma;

  c_1 = params.c_1;
  cx_1 = params.cx_1;
  phix_1 = params.phix_1;
  eta_1 = params.eta_1;

  % At x = 1
  alpha_3_N = -alpha*exp(-alpha*n*eta_1) - (1-alpha)*exp((1-alpha)*n*eta_1);
  k_c_approx = sqrt(1./(alpha_3_N*gamma).*(-alpha_3_N*n.*phix_1 + exp(-alpha*n*eta_1).*cx_1./(c_1-rho_s)));

end

function [k_max_approx, omega_max_approx] = get_max_approx(params, k_c_approx)

  [k_max_approx, omega_max_approx_temp] = fminbnd(@(ks_i) -get_omega_approx(params, ks_i), 0, k_c_approx);
  omega_max_approx = -omega_max_approx_temp;

end

function omega_approx = get_omega_approx(params, ks_i)

  z_p = params.z_p;
  z = params.z;
  D_p = params.D_p;
  D = params.D;
  rho_s = params.rho_s;
  n = params.n;
  alpha = params.alpha;
  gamma = params.gamma;
  beta_m = params.beta_m;
  beta_v = params.beta_v;

  c_1 = params.c_1;
  cx_1 = params.cx_1;
  phix_1 = params.phix_1;
  j_0_1 = params.j_0_1;
  eta_1 = params.eta_1;
  ct_1 = params.ct_1;

  alpha_1 = D - D_p;
  alpha_2 = z_p*D_p - z*D;

  c_1_i = c_1;
  cx_1_i = cx_1;
  phix_1_i = phix_1;
  j_0_1_i = j_0_1;
  eta_1_i = eta_1;
  ct_1_i = ct_1;

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

  omega_approx = omega_approx_i;

end
