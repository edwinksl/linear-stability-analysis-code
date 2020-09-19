function varargout = get_max_approx(case_flag, save_flag)

  elapsed_time = tic;

  bc_index = 1;

  if strcmp(case_flag, 'ca')

    base_dir = 'datasets/convergence-analysis/base/';
    lsa_dir = 'datasets/convergence-analysis/lsa/';

    rho_ss = [-0.05, 0, 0.05];
    Das = 1;
    J_applieds = 1.5;
    Ns = [251, 501, 1001, 2001, 4001];

  elseif strcmp(case_flag, 'a')

    base_dir = 'datasets/analysis/base/';
    lsa_dir = 'datasets/analysis/lsa/';

    rho_ss = [-0.05, 0, 0.05];
    Das = [0.1, 1, 10];
    J_applieds = [0.5, 1, 1.5];
    Ns = 1001;

  elseif strcmp(case_flag, 'ps')

    base_dir = 'datasets/parameter-sweeps/base/';
    lsa_dir = 'datasets/parameter-sweeps/lsa/';

    rho_ss = [-1, -0.75, -0.5, -0.25, -0.05];
    Das = 1;
    J_applieds = 1.5;
    Ns = 1001;

  end

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_N = length(Ns);

  k_max_approx = cell(N_rho_s, N_Da, N_J_applied, N_N);
  omega_max_approx = cell(N_rho_s, N_Da, N_J_applied, N_N);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        for l = 1:N_N

          N = Ns(l);

          str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
          base = lsa_base(base_dir, str, bc_index);

          N_t = length(base.t);
          t_indices = 1:N_t;
          k_c_0 = get_k_c_approx(base, t_indices);

          k_max = NaN(N_t, 1);
          omega_max = NaN(N_t, 1);

          parfor m = 1:N_t

            [k_max_m, omega_max_temp_m] = fminbnd(@(ks) -obj_func(base, str, bc_index, m, ks), 0, k_c_0(m));

            k_max(m) = k_max_m;
            omega_max(m) = -omega_max_temp_m;

          end

          k_max_approx{i,j,k,l} = k_max;
          omega_max_approx{i,j,k,l} = omega_max;

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

  omega_approx = get_omega_approx(base, str, bc_index, t_indices, {ks});
  y = omega_approx{1};

end
