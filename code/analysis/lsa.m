function varargout = lsa(base, str, bc_index, t_indices, ks, solver_index, save_flag, lsa_dir, varargin)

  % For brevity, drop subscript 0 that denotes base state

  % Use varargin to pass in flag for pulse charging

  % Assume bc_index = 1 throughout

  % t_indices = indices of t

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
  % beta_D = base.beta_D;
  beta_v = base.beta_v;
  Da = base.Da;
  t_s = base.t_s;
  J_applied = base.J_applied;
  xi_p = base.xi_p;
  dx = base.dx;

  x = base.x;
  N = length(x);
  t = base.t;
  t_ratio = t/t_s;

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
    assert(parameters(5) == delta_t_on/t_s)
    assert(parameters(6) == gamma_dc)
  else
    t_correct = round(get_t_ratio(parameters), 6);
    t_base = round(t_ratio, 6);
    assert(isequal(t_correct, t_base))
  end

  % fprintf('z_p = %d\n', z_p)
  % fprintf('z = %d\n', z)
  % fprintf('D_p = %g\n', D_p)
  % fprintf('D = %g\n', D)
  fprintf('rho_s = %g\n', rho_s)
  % fprintf('n = %d\n', n)
  % fprintf('alpha = %g\n', alpha)
  % fprintf('gamma = %g\n', gamma)
  % fprintf('beta_m = %g\n', beta_m)
  % fprintf('beta_D = %g\n', beta_D)
  % fprintf('beta_v = %g\n', beta_v)
  fprintf('Da = %g\n', Da)
  % fprintf('t_s = %g\n', t_s)
  fprintf('J_applied = %g\n', J_applied)
  % fprintf('xi_p = %g\n', xi_p)
  % fprintf('dx = %g\n', dx)
  fprintf('N = %d\n', N)
  if pulse_flag == 1
    fprintf('delta_t_on/t_s = %g\n', delta_t_on/t_s)
    fprintf('gamma_dc = %g\n', gamma_dc)
  end
  fprintf('\n')

  c = base.c;
  cx = base.cx;
  phix = base.phix;
  phixx = base.phixx;

  j_0_0 = base.j_0_0;
  j_0_1 = base.j_0_1;
  eta_0 = base.eta_0;
  eta_1 = base.eta_1;
  ct_0 = base.ct_0;
  ct_1 = base.ct_1;

  % Parameters
  alpha_1 = D - D_p;
  alpha_2 = z_p*D_p - z*D;

  N_tot = 2*N + 2;

  N_t = length(t_indices);
  omega_mu = cell(N_t, 1);
  % v_mu = cell(N_t, 1);
  h_1_a_mu = cell(N_t, 1);
  h_1_c_mu = cell(N_t, 1);

  % Construct Z, which is independent of t and k
  Z_m_2 = zeros(N_tot, 1);
  Z_m_2(1) = 1;

  Z_p_2 = zeros(N_tot, 1);
  Z_p_2(end) = -1;

  Z_0_temp = repmat([1;0], [N-2,1]);

  Z_0 = zeros(N_tot, 1);
  Z_0(4:end-3) = Z_0_temp;

  if solver_index == 1

    Z_0(1) = 1;
    Z_0(end) = -1;

  end

  Z = spdiags([Z_m_2, Z_0, Z_p_2], [-2, 0, 2], N_tot, N_tot);

  if solver_index == 1

    Z = full(Z);
    perm_indices = NaN;

  elseif solver_index == 2

    perm_indices = NaN(N_tot, 1);
    perm_indices(1) = 3;
    for i = 2:N-1
      perm_indices(i) = 2*i;
    end
    perm_indices(N) = 2*N;
    perm_indices(N+1) = 1;
    perm_indices(N+2) = 2;
    for i = N+3:2*N
      perm_indices(i) = 2*(i-(N+3)) + 5;
    end
    perm_indices(2*N+1) = 2*N + 1;
    perm_indices(2*N+2) = 2*N + 2;

    Z = Z(perm_indices, :);

  end

  parfor i = 1:N_t

    ks_i = ks{i};
    N_k_i = length(ks_i);
    omega_mu_i = NaN(N_k_i, 1);
    % v_mu_i = NaN(N_tot, N_k_i);
    h_1_a_mu_i = NaN(N_k_i, 1);
    h_1_c_mu_i = NaN(N_k_i, 1);

    t_temp = t;
    c_temp = c;
    cx_temp = cx;
    phix_temp = phix;
    phixx_temp = phixx;

    j_0_0_temp = j_0_0;
    j_0_1_temp = j_0_1;
    eta_0_temp = eta_0;
    eta_1_temp = eta_1;
    ct_0_temp = ct_0;
    ct_1_temp = ct_1;

    t_index = t_indices(i);

    t_i = t_temp(t_index);

    if t_i == 0

      beta_1 = 1 + (rho_s+abs(rho_s))/2;

      % For all x
      c_i = beta_1*ones(N, 1);
      cx_i = zeros(N, 1);
      phix_i = zeros(N, 1);
      phixx_i = zeros(N, 1);

      j_0_0_i = Da*n*(xi_p*(beta_1-rho_s))^(1-alpha);
      j_0_1_i = Da*n*(xi_p*(beta_1-rho_s))^(1-alpha);
      eta_0_i = 0;
      eta_1_i = 0;
      ct_0_i = 0;
      ct_1_i = 0;

    else

      % For all x
      c_i = c_temp(:, t_index);
      cx_i = cx_temp(:, t_index);
      phix_i = phix_temp(:, t_index);
      phixx_i = phixx_temp(:, t_index);

      j_0_0_i = j_0_0_temp(t_index);
      j_0_1_i = j_0_1_temp(t_index);
      eta_0_i = eta_0_temp(t_index);
      eta_1_i = eta_1_temp(t_index);
      ct_0_i = ct_0_temp(t_index);
      ct_1_i = ct_1_temp(t_index);

    end

    % At x = 0
    c_0 = c_i(1);  % c at x = 0
    cx_0 = cx_i(1);  % cx at x = 0
    phix_0 = phix_i(1);  % phix at x = 0

    alpha_3_1 = -alpha*exp(-alpha*n*eta_0_i) - (1-alpha)*exp((1-alpha)*n*eta_0_i);  % alpha_3 at x = 0
    % alpha_4_1 = exp(-alpha*n*eta_0_i) - exp((1-alpha)*n*eta_0_i);  % alpha_4 at x = 0

    B_temp = alpha_2*c_0 - z_p*D_p*rho_s;
    D_temp = beta_v*j_0_0_i;

    % At x = 1
    c_1 = c_i(end);  % c at x = 1
    cx_1 = cx_i(end);  % cx at x = 1
    phix_1 = phix_i(end);  % phix at x = 1

    alpha_3_N = -alpha*exp(-alpha*n*eta_1_i) - (1-alpha)*exp((1-alpha)*n*eta_1_i);  % alpha_3 at x = 1
    % alpha_4_N = exp(-alpha*n*eta_1_i) - exp((1-alpha)*n*eta_1_i);  % alpha_4 at x = 1

    E_temp = alpha_2*c_1 - z_p*D_p*rho_s;
    G_temp = beta_v*j_0_1_i;

    % Interior values
    c_int = c_i(2:end-1);
    cx_int = cx_i(2:end-1);
    phix_int = phix_i(2:end-1);
    phixx_int = phixx_i(2:end-1);

    % Construct t-dependent and k-independent entries of Y
    M1 = D*(1/dx^2 - z*phix_int/(2*dx));
    M3 = D*(1/dx^2 + z*phix_int/(2*dx));
    M4 = z*D*(-cx_int/(2*dx) + c_int/dx^2);
    M6 = z*D*(cx_int/(2*dx) + c_int/dx^2);

    A1 = alpha_1*1/dx^2 + alpha_2*phix_int/(2*dx);
    A3 = alpha_1*1/dx^2 - alpha_2*phix_int/(2*dx);
    A4 = z_p*D_p*rho_s*1/dx^2 - alpha_2*(-cx_int/(2*dx)+c_int/dx^2);
    A6 = z_p*D_p*rho_s*1/dx^2 - alpha_2*(cx_int/(2*dx)+c_int/dx^2);

    B1 = beta_m*(alpha_1*3/(2*dx) + alpha_2*phix_0);
    B2 = beta_m*(-alpha_1*2/dx);
    B3 = beta_m*(alpha_1*1/(2*dx));
    B4 = beta_m*B_temp*(-3/(2*dx));
    B5 = beta_m*B_temp*(2/dx);
    B6 = beta_m*B_temp*(-1/(2*dx));

    C1 = ct_0_i;
    C2 = -D*(3/(2*dx) - z*phix_0);
    C3 = -D*(-2/dx);
    C4 = -D*(1/(2*dx));
    C5 = z*D*c_0*(-3/(2*dx));
    C6 = z*D*c_0*(2/dx);
    C7 = z*D*c_0*(-1/(2*dx));

    D2 = D_temp*(exp(-alpha*n*eta_0_i)/(c_0-rho_s));
    D3 = D_temp*(-alpha_3_1*n);

    E1 = beta_m*(alpha_1*1/(2*dx));
    E2 = beta_m*(-alpha_1*2/dx);
    E3 = beta_m*(alpha_1*3/(2*dx) - alpha_2*phix_1);
    E4 = beta_m*E_temp*(-1/(2*dx));
    E5 = beta_m*E_temp*(2/dx);
    E6 = beta_m*E_temp*(-3/(2*dx));

    F1 = -ct_1_i;
    F2 = -D*(1/(2*dx));
    F3 = -D*(-2/dx);
    F4 = -D*(3/(2*dx) + z*phix_1);
    F5 = z*D*c_1*(-1/(2*dx));
    F6 = z*D*c_1*(2/dx);
    F7 = z*D*c_1*(-3/(2*dx));

    G2 = G_temp*(exp(-alpha*n*eta_1_i)/(c_1-rho_s));
    G3 = G_temp*(-alpha_3_N*n);

    Y_m_5 = zeros(N_tot, 1);
    Y_m_5(end-6) = F2;

    Y_m_4 = zeros(N_tot, 1);
    Y_m_4(end-6) = E1;
    Y_m_4(end-5) = F5;

    Y_m_3 = zeros(N_tot, 1);
    Y_m_3(2:2:2*N-4) = A1;
    Y_m_3(end-5) = E4;
    Y_m_3(end-4) = F3;

    Y_m_2 = zeros(N_tot, 1);
    Y_m_2(2:2:2*N-4) = M1;
    Y_m_2(3:2:2*N-3) = A4;
    Y_m_2(end-4) = E2;
    Y_m_2(end-3) = F6;
    Y_m_2(end-2) = G2;

    Y_p_2 = zeros(N_tot, 1);
    Y_p_2(3) = D3;
    Y_p_2(4) = C3;
    Y_p_2(5) = B5;
    Y_p_2(6:2:end-2) = M3;
    Y_p_2(7:2:end-1) = A6;

    Y_p_3 = zeros(N_tot, 1);
    Y_p_3(5) = C6;
    Y_p_3(6) = B3;
    Y_p_3(7:2:end-1) = M6;

    Y_p_4 = zeros(N_tot, 1);
    Y_p_4(6) = C4;
    Y_p_4(7) = B6;

    Y_p_5 = zeros(N_tot, 1);
    Y_p_5(7) = C7;

    % if solver_index == 1

      % Hack to fix "[u]ndefined function or variable" error messages
      G7_p = NaN;
      D2_p = NaN;
      Y_m_6 = NaN;
      Y_p_6 = NaN;
      % perm_indices = NaN;

    % elseif solver_index == 2
    if solver_index == 2

      D2_p = D2 - B1;
      D3_p = -B2;
      D4_p = -B3;
      D5_p = D3 - B4;
      D6_p = -B5;
      D7_p = -B6;

      G2_p = -E1;
      G3_p = -E2;
      G4_p = G2 - E3;
      G5_p = -E4;
      G6_p = -E5;
      G7_p = G3 - E6;

      Y_m_6 = zeros(N_tot, 1);
      Y_m_6(end-6) = G2_p;

      Y_m_5(end-5) = G5_p;

      Y_m_4(end-4) = G3_p;

      Y_m_3(end-3) = G6_p;

      Y_m_2(end-2) = G4_p;

      Y_p_2(3) = D5_p;

      Y_p_3(4) = D3_p;

      Y_p_4(5) = D6_p;

      Y_p_5(6) = D4_p;

      Y_p_6 = zeros(N_tot, 1);
      Y_p_6(7) = D7_p;

    end

    for k_index = 1:N_k_i

      k = ks_i(k_index);

      % Construct t-dependent and k-dependent entries of Y
      M2 = D*(-2/dx^2 - k^2 + z*phixx_int);
      M5 = z*D*(-2*c_int/dx^2 - k^2*c_int);
      A2 = alpha_1*(-2/dx^2-k^2) - alpha_2*phixx_int;
      A5 = z_p*D_p*rho_s*(-2/dx^2-k^2) - alpha_2*(-2*c_int/dx^2-k^2*c_int);
      D1 = D_temp*(alpha_3_1*n*(-phix_0+gamma*k^2/n) + exp(-alpha*n*eta_0_i)*cx_0/(c_0-rho_s));
      G1 = G_temp*(alpha_3_N*n*(-phix_1-gamma*k^2/n) + exp(-alpha*n*eta_1_i)*cx_1/(c_1-rho_s));

      Y_m_1 = zeros(N_tot, 1);
      Y_m_1(1) = C1;
      Y_m_1(2) = B1;
      Y_m_1(3:2:2*N-3) = M4;
      Y_m_1(4:2:2*N-2) = A2;
      Y_m_1(end-3) = E5;
      Y_m_1(end-2) = F4;
      Y_m_1(end-1) = G3;

      Y_0 = zeros(N_tot, 1);
      Y_0(1) = D1;
      Y_0(2) = C2;
      Y_0(3) = B4;
      Y_0(4:2:2*N-2) = M2;
      Y_0(5:2:2*N-1) = A5;
      Y_0(end-2) = E3;
      Y_0(end-1) = F7;
      Y_0(end) = G1;

      Y_p_1 = zeros(N_tot, 1);
      Y_p_1(2) = D2;
      Y_p_1(3) = C5;
      Y_p_1(4) = B2;
      Y_p_1(5:2:end-3) = M5;
      Y_p_1(6:2:end-2) = A3;
      Y_p_1(end-1) = E6;
      Y_p_1(end) = F1;

      if solver_index == 2

        D1_p = D1;
        G1_p = G1;

        Y_m_1(end-1) = G7_p;

        Y_0(1) = D1_p;
        Y_0(end) = G1_p;

        Y_p_1(2) = D2_p;

      end

      if solver_index == 1

        Y = spdiags([Y_m_5, Y_m_4, Y_m_3, Y_m_2, Y_m_1, Y_0, Y_p_1, Y_p_2, Y_p_3, Y_p_4, Y_p_5], -5:5, N_tot, N_tot);

        [V, omega] = eig(full(Y), Z, 'vector');
        % omega = eig(full(Y), Z);

        % Check numbers of finite, infinite and NaN eigenvalues
        N_finite = sum(isfinite(omega));
        N_inf = sum(isinf(omega));
        N_NaN = sum(isnan(omega));

        omega = omega(isfinite(omega));

        % Check numbers of real and complex finite eigenvalues
        N_real = sum(isreal(omega));
        N_complex = sum(~isreal(omega));

        % omega

        [~, omega_mu_index] = max(real(omega));  % mu = most unstable
        omega_mu_ki = omega(omega_mu_index);
        v_mu_ki = V(:, omega_mu_index);
        h_1_a = v_mu_ki(1);
        h_1_c = v_mu_ki(end);
        omega_mu_i(k_index) = omega_mu_ki;
        % v_mu_i(:, k_index) = v_mu_ki;
        h_1_a_mu_i(k_index) = h_1_a;
        h_1_c_mu_i(k_index) = h_1_c;

        t_is = t_i/t_s;
        % h_ratio = h_1_c/h_1_a;
        % fprintf('t/t_s = %g, k = %g, rho_s = %g, Da = %g, J_applied = %g, N = %d, omega_mu = %g+%gi, h_1_a = %g+%gi, h_1_c = %g+%gi, h_ratio = %g+%gi, N_finite = %d, N_real = %d, N_complex = %d, N_inf = %d, N_NaN = %d, solver_index = %d\n', t_is, k, rho_s, Da, J_applied, N, real(omega_mu_ki), imag(omega_mu_ki), real(h_1_a), imag(h_1_a), real(h_1_c), imag(h_1_c), real(h_ratio), imag(h_ratio), N_finite, N_real, N_complex, N_inf, N_NaN, solver_index)
        fprintf('t/t_s = %g, k = %g, rho_s = %g, Da = %g, J_applied = %g, N = %d, omega_mu = %g+%gi, N_finite = %d, N_real = %d, N_complex = %d, N_inf = %d, N_NaN = %d, solver_index = %d\n', t_is, k, rho_s, Da, J_applied, N, real(omega_mu_ki), imag(omega_mu_ki), N_finite, N_real, N_complex, N_inf, N_NaN, solver_index)

      elseif solver_index == 2

        Y = spdiags([Y_m_6, Y_m_5, Y_m_4, Y_m_3, Y_m_2, Y_m_1, Y_0, Y_p_1, Y_p_2, Y_p_3, Y_p_4, Y_p_5, Y_p_6], -6:6, N_tot, N_tot);

        Y = Y(perm_indices, :);

        Z_hat = Z(1:N, :);
        sigma = -1e2;
        E = -sigma*eye(N+2, N+2);
        H = Y(1:N, :);
        P = Y(N+1:end, :);
        Y_bar = [H; E*P];
        Z_bar = [Z_hat; -P];

        t_is = t_i/t_s;
        [N_dim, tol, V, omega, conv_flag] = eigs_wrapper(N, t_is, k, rho_s, Da, J_applied, Y_bar, Z_bar);
        omega = diag(omega);

        % Check numbers of finite, infinite and NaN eigenvalues
        N_finite = sum(isfinite(omega));
        N_inf = sum(isinf(omega));
        N_NaN = sum(isnan(omega));

        % Check numbers of real and complex finite eigenvalues
        N_real = sum(isreal(omega));
        N_complex = sum(~isreal(omega));

        % omega

        [~, omega_mu_index] = max(real(omega));  % mu = most unstable
        omega_mu_ki = omega(omega_mu_index);
        v_mu_ki = V(:, omega_mu_index);
        h_1_a = v_mu_ki(1);
        h_1_c = v_mu_ki(end);
        omega_mu_i(k_index) = omega_mu_ki;
        % v_mu_i(:, k_index) = v_mu_ki;
        h_1_a_mu_i(k_index) = h_1_a;
        h_1_c_mu_i(k_index) = h_1_c;

        % h_ratio = h_1_c/h_1_a;
        % fprintf('t/t_s = %g, k = %g, rho_s = %g, Da = %g, J_applied = %g, N = %d, omega_mu = %g+%gi, h_1_a = %g+%gi, h_1_c = %g+%gi, h_ratio = %g+%gi, N_finite = %d, N_real = %d, N_complex = %d, N_inf = %d, N_NaN = %d, N_dim = %d, tol = %g, solver_index = %d, conv_flag = %d\n', t_is, k, rho_s, Da, J_applied, N, real(omega_mu_ki), imag(omega_mu_ki), real(h_1_a), imag(h_1_a), real(h_1_c), imag(h_1_c), real(h_ratio), imag(h_ratio), N_finite, N_real, N_complex, N_inf, N_NaN, N_dim, tol, solver_index, conv_flag)
        fprintf('t/t_s = %g, k = %g, rho_s = %g, Da = %g, J_applied = %g, N = %d, omega_mu = %g+%gi, N_finite = %d, N_real = %d, N_complex = %d, N_inf = %d, N_NaN = %d, N_dim = %d, tol = %g, solver_index = %d, conv_flag = %d\n', t_is, k, rho_s, Da, J_applied, N, real(omega_mu_ki), imag(omega_mu_ki), N_finite, N_real, N_complex, N_inf, N_NaN, N_dim, tol, solver_index, conv_flag)

      end

    end

    fprintf('\n')

    omega_mu{i} = omega_mu_i;
    % v_mu{i} = v_mu_i;
    h_1_a_mu{i} = h_1_a_mu_i;
    h_1_c_mu{i} = h_1_c_mu_i;

  end

  varargout{1} = omega_mu;
  varargout{2} = base;
  % varargout{3} = v_mu;
  varargout{3} = h_1_a_mu;
  varargout{4} = h_1_c_mu;

  if save_flag == 1
    filename = ['lsa', '_', num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];

    if bc_index == 1
      subdir = 'applied-current/';
    end

    lsa_full_dir = [lsa_dir, subdir];

    % save([lsa_full_dir, filename, '.mat'], 'omega_mu', 'base', 'v_mu', 'bc_index', 't_indices', 'ks', 'solver_index')
    save([lsa_full_dir, filename, '.mat'], 'omega_mu', 'base', 'h_1_a_mu', 'h_1_c_mu', 'bc_index', 't_indices', 'ks', 'solver_index')
  end

end

function [N_dim, tol] = get_parameters(N, t_is, k, rho_s, Da, J_applied)

  if N == 251
    N_dim = 25;
    tol = 2.5e-9;
  elseif N == 501
    N_dim = 35;
    tol = 2.5e-9;
  elseif N == 1001 || N == 2001 || N == 4001
    if N == 1001
      N_dim = 65;
    elseif N == 2001
      N_dim = 65;
    elseif N == 4001
      N_dim = 70;
    end
    if t_is == 0
      tol = 1e-12;
    else
      tol = 5e-9;
    end
  end

end

function t_ratio = get_t_ratio(parameters)

  rho_s = parameters(1);
  Da = parameters(2);
  J_applied = parameters(3);
  N = parameters(4);

  t_all = [0:0.2:0.8, 0.85:0.05:0.95, 1:0.2:2];

  if rho_s < 0 || ((rho_s == 0 || rho_s == 0.05) && J_applied == 0.5)

    t_ratio = t_all;

  elseif rho_s == 0 && (J_applied == 1.5 || J_applied == 1)

    t_ratio = t_all(1:8);

  elseif rho_s == 0.05 && (J_applied == 1.5 || J_applied == 1)

    t_ratio = t_all(1:6);

  end

end

function [V, omega] = expm_wrapper(Y_bar, Z_bar)

  addpath('expLReigs')

  % Hack to suppress fprintf statements; see https://stackoverflow.com/q/3029636/486919
  [T, V, omega, eig_res] = evalc('expm_LR_eigs(2, Y_bar, Z_bar, 1)');

end

function [N_dim, tol, V, omega, conv_flag] = eigs_wrapper(N, t_is, k, rho_s, Da, J_applied, Y_bar, Z_bar)

  if k < 2
    N_dim = NaN;
    tol = NaN;
    [V, omega] = expm_wrapper(Y_bar, Z_bar);
    conv_flag = NaN;
  else
    N_omega = 1;
    [N_dim, tol] = get_parameters(N, t_is, k, rho_s, Da, J_applied);
    [V, omega, conv_flag] = eigs(Y_bar, Z_bar, N_omega, 'largestreal', 'SubspaceDimension', N_dim, 'Tolerance', tol);
    if conv_flag == 1
      [V, omega] = expm_wrapper(Y_bar, Z_bar);
    end
  end

end
