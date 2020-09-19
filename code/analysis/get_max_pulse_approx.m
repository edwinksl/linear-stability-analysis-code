function varargout = get_max_pulse_approx(sign_flag, save_flag)

  elapsed_time = tic;

  bc_index = 1;

  if strcmp(sign_flag, 'm')

    base_dir = '/home/edwinksl/Dropbox/git/linear-stability-analysis/code/analysis/datasets/pulse/base/negative/';
    lsa_dir = 'datasets/pulse/lsa/negative/';
    rho_ss = -0.05;
    Das = 1;
    J_applieds = 1.5;  % average J_applied, not peak J_applied
    N = 1001;
    delta_t_ons = 1;  % delta_t_on/t_s, not just delta_t_on
    gamma_dcs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
    tol = 5e-1;

  elseif strcmp(sign_flag, 'z')

    base_dir = '/home/edwinksl/Dropbox/git/linear-stability-analysis/code/analysis/datasets/pulse/base/zero/';
    lsa_dir = 'datasets/pulse/lsa/zero/';
    rho_ss = 0;
    Das = 1;
    J_applieds = 1;  % average J_applied, not peak J_applied
    N = 1001;
    delta_t_ons = 0.0125;  % delta_t_on/t_s, not just delta_t_on
    gamma_dcs = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
    tol = 5e-1;

  end

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_delta_t_on = length(delta_t_ons);
  N_gamma_dc = length(gamma_dcs);

  k_max_approx = cell(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);
  omega_max_approx = cell(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        for l = 1:N_delta_t_on

          delta_t_on = delta_t_ons(l);

          for m = 1:N_gamma_dc

            gamma_dc = gamma_dcs(m);

            str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N), '_', num2str(delta_t_on), '_', num2str(gamma_dc)];
            base = lsa_base(base_dir, str, bc_index);

            N_t = length(base.t);
            t_indices = 1:N_t;
            k_c_0 = get_k_c_approx(base, t_indices);

            k_c_0 = remove_imag(k_c_0, tol);

            k_max = NaN(N_t, 1);
            omega_max = NaN(N_t, 1);

            parfor n = 1:N_t

              [k_max_n, omega_max_temp_n] = fminbnd(@(ks) -obj_func(base, str, bc_index, n, ks), 0, k_c_0(n));

              k_max(n) = k_max_n;
              omega_max(n) = -omega_max_temp_n;

            end

            k_max_approx{i,j,k,l,m} = k_max;
            omega_max_approx{i,j,k,l,m} = omega_max;

          end

        end

      end

    end

  end

  varargout{1} = k_max_approx;
  varargout{2} = omega_max_approx;

  if save_flag == 1
    filename = 'max_approx';

    if bc_index == 1
      subdir = 'applied-current/';
    end

    lsa_full_dir = [lsa_dir, subdir];

    save([lsa_full_dir, filename, '.mat'], 'k_max_approx', 'omega_max_approx')
  end

  total_elapsed_time = toc(elapsed_time);
  fprintf('Total elapsed time = %g s\n', total_elapsed_time)

end

function y = obj_func(base, str, bc_index, t_indices, ks)

  omega_approx = get_omega_approx(base, str, bc_index, t_indices, {ks}, 1);
  y = omega_approx{1};

end

function x = remove_imag(x, tol)

  indices = imag(x) < tol;
  x_0 = x(indices);
  x(indices) = real(x_0);

end
