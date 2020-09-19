function varargout = get_c_max_pulse_num(save_flag)

  elapsed_time = tic;

  bc_index = 1;
  solver_index = 2;
  base_dir = 'datasets/pulse/base/';
  lsa_dir = 'datasets/pulse/lsa/';

  rho_ss = -0.05;
  Das = 1;
  J_applieds = 1.5;  % average J_applied, not peak J_applied
  N = 1001;
  delta_t_ons = 1;  % delta_t_on/delta_t_s, not just delta_t_on
  gamma_dcs = 0.1;

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_delta_t_on = length(delta_t_ons);
  N_gamma_dc = length(gamma_dcs);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;

  k_c_num = cell(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);
  k_max_num = cell(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);
  omega_max_num = cell(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);

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

            t = base.t;
            N_t = length(t);
            t_indices = 1:N_t;
            k_c_0 = get_k_c_approx(base, t_indices);

            k_c = NaN(N_t, 1);
            options_k_c = optimset('TolX', 1e-6);

            k_max = NaN(N_t, 1);
            omega_max = NaN(N_t, 1);

            k_max_approx_ijklm = k_max_approx{i,j,k,l,m};

            parfor n = 1:N_t

              t_temp = t;

              t_index = t_indices(n);

              t_n = t_temp(t_index);

              if t_n == 0
                % k_c(n) = fzero(@(ks) obj_func(base, str, bc_index, n, ks, solver_index, save_flag, lsa_dir), [0, k_c_0(n)], options_k_c);
                k_c(n) = 0;
              else
                k_c(n) = fzero(@(ks) obj_func(base, str, bc_index, n, ks, solver_index, save_flag, lsa_dir), k_c_0(n), options_k_c);
              end

              if t_n == 0
                k_max_n = 0;
                omega_max_temp_n = -0;
              else
                [k_max_n, omega_max_temp_n] = fminbnd(@(ks) -obj_func(base, str, bc_index, n, ks, solver_index, save_flag, lsa_dir), 0, k_max_approx_ijklm(n));
              end

              k_max(n) = k_max_n;
              omega_max(n) = -omega_max_temp_n;

            end

            k_c_num{i,j,k,l,m} = k_c;
            k_max_num{i,j,k,l,m} = k_max;
            omega_max_num{i,j,k,l,m} = omega_max;

          end

        end

      end

    end

  end

  varargout{1} = k_c_num;
  varargout{2} = k_max_num;
  varargout{3} = omega_max_num;

  if save_flag == 1
    filename = 'c_max_num';

    if bc_index == 1
      subdir = 'applied-current/';
    end

    lsa_full_dir = [lsa_dir, subdir];

    save([lsa_full_dir, filename, '.mat'], 'k_c_num', 'k_max_num', 'omega_max_num')
  end

  total_elapsed_time = toc(elapsed_time);
  fprintf('Total elapsed time = %g s\n', total_elapsed_time)

end

function y = obj_func(base, str, bc_index, t_indices, ks, solver_index, ~, lsa_dir)

  omega_mu = lsa(base, str, bc_index, t_indices, {ks}, solver_index, 0, lsa_dir, 1);
  y = omega_mu{1};

end
