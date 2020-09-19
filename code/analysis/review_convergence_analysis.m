function review_convergence_analysis(save_flag)

  close all

  bc_index = 1;
  base_dir = 'datasets/convergence-analysis/base/';
  lsa_dir = 'datasets/convergence-analysis/lsa/';
  graphics_dir = '../../review-graphics/';

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  rho_ss = 0;
  Das = 1;
  J_applieds = 1.5;
  Ns = 4001;

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_N = length(Ns);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;
  omega_max_approx = max_approx_data.omega_max_approx;

  c_max_num_data = load([lsa_full_dir, 'c_max_num', '.mat']);
  k_c_num = c_max_num_data.k_c_num;
  k_max_num = c_max_num_data.k_max_num;
  omega_max_num = c_max_num_data.omega_max_num;

  k_c_approx_t = cell(N_rho_s, N_Da, N_J_applied, N_N);
  k_max_approx_t = cell(N_rho_s, N_Da, N_J_applied, N_N);
  omega_max_approx_t = cell(N_rho_s, N_Da, N_J_applied, N_N);
  k_c_num_t = cell(N_rho_s, N_Da, N_J_applied, N_N);
  k_max_num_t = cell(N_rho_s, N_Da, N_J_applied, N_N);
  omega_max_num_t = cell(N_rho_s, N_Da, N_J_applied, N_N);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        ks_combined = cell(N_N, 1);
        log10_ks_combined = cell(N_N, 1);
        omega_mu_combined = cell(N_N, 1);
        omega_approx_combined = cell(N_N, 1);
        legends_combined = cell(2*N_N, 1);

        omega_mu_k_combined = cell(N_N, 1);
        omega_mu_kk_combined = cell(N_N, 1);
        omega_approx_k_combined = cell(N_N, 1);
        omega_approx_kk_combined = cell(N_N, 1);

        for l = 1:N_N

          N = Ns(l);

          str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
          filename = ['lsa', '_', str];
          lsa_data = load([lsa_full_dir, filename, '.mat']);

          base = lsa_data.base;
          ks = lsa_data.ks;
          omega_mu = lsa_data.omega_mu;
          t_index = lsa_data.t_index;

          t = base.t;
          t_s = base.t_s;
          rho_s = base.rho_s;
          Da = base.Da;
          J_applied = base.J_applied;

          t_rs = t/t_s;  % rs = rescaled
          t_rs_num = t_rs(t_index);
          log10_ks = cellfun(@log10, ks, 'UniformOutput', false);
          omega_approx = get_omega_approx(base, str, bc_index, t_index, ks);

          k_c_approx_t{i,j,k,l} = get_k_c_approx(base, t_index);
          k_max_approx_t{i,j,k,l} = k_max_approx{i,j,k,l}(t_index);
          omega_max_approx_t{i,j,k,l} = omega_max_approx{i,j,k,l}(t_index);
          k_c_num_t{i,j,k,l} = k_c_num{i,j,k,l}(t_index);
          k_max_num_t{i,j,k,l} = k_max_num{i,j,k,l}(t_index);
          omega_max_num_t{i,j,k,l} = omega_max_num{i,j,k,l}(t_index);

          % Index into cells
          ks = ks{1};
          log10_ks = log10_ks{1};
          omega_mu = omega_mu{1};
          omega_approx = omega_approx{1};

          ks_combined{l} = ks;
          log10_ks_combined{l} = log10_ks;
          omega_mu_combined{l} = omega_mu;
          omega_approx_combined{l} = omega_approx;

          legends_combined{l} = ['$$N$$ = ', num2str(N), ' (num.)'];
          legends_combined{l+N_N} = ['$$N$$ = ', num2str(N), ' (approx.)'];

          omega_mu_k = gradient(omega_mu, ks);
          omega_mu_kk = gradient(omega_mu_k, ks);
          omega_approx_k = gradient(omega_approx, ks);
          omega_approx_kk = gradient(omega_approx_k, ks);

          omega_mu_k_combined{l} = omega_mu_k;
          omega_mu_kk_combined{l} = omega_mu_kk;
          omega_approx_k_combined{l} = omega_approx_k;
          omega_approx_kk_combined{l} = omega_approx_kk;

        end

        title_combined_str = ['$$t/t_\textnormal{s}$$ = ', num2str(t_rs_num, '%.3g'), ', ', '$$\rho_\textnormal{s}$$ = ', num2str(rho_s, '%.3g'), ', ', '$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$J_\textnormal{a}$$ = ', num2str(J_applied, '%.3g')];  % rename J_applied to J_a in title; use t_rs_num from last iteration of for loop

        save_str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied)];

        if save_flag == 1
          fig = figure('Visible', 'off');
        else
          fig = figure;
        end

        fig.Renderer = 'painters';

        ks_combined_linear = horzcat(ks_combined{:});
        omega_mu_combined_linear = vertcat(omega_mu_combined{:});
        omega_approx_combined_linear = horzcat(omega_approx_combined{:});
        omega_mu_real = real(omega_mu_combined_linear);
        omega_approx_real = real(omega_approx_combined_linear);

        k_min = min(ks_combined_linear);
        k_max = max(ks_combined_linear);
        omega_min = min(min(omega_mu_real), min(omega_approx_real));
        omega_max = max(max(omega_mu_real), max(omega_approx_real));

        patch([k_min, k_max, k_max, k_min], [0, 0, omega_max, omega_max], 'red', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')
        patch([k_min, k_max, k_max, k_min], [0, 0, omega_min, omega_min], 'blue', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')

        ax = gca;

        ax.ColorOrderIndex = 5;

        for l = 1:N_N
          hold on
          plot(ks_combined{l}, omega_mu_combined{l}, 'LineWidth', 2)
          hold off
        end

        ax.ColorOrderIndex = 5;

        for l = 1:N_N
          hold on
          plot(ks_combined{l}, omega_approx_combined{l}, '--', 'LineWidth', 2)
          hold off
        end

        xlabel('$$k$$', 'Interpreter', 'latex')
        ylabel('$$\Re\left(\omega\right)$$', 'Interpreter', 'latex')
        lgd = legend(legends_combined, 'Location', 'south', 'Interpreter', 'latex');
        lgd.NumColumns = 2;
        title(title_combined_str, 'Interpreter', 'latex')
        axis([k_min, k_max, omega_min, omega_max])
        ax.FontSize = 12;
        ax.TickDir = 'out';  % tick marks directed inward are not visible in EPS files for unknown reasons

        if save_flag == 1
          filename = ['review_ca', '_', save_str];
          saveas(gcf, [graphics_dir, filename, '.eps'], 'epsc')
        end

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;

        ax.ColorOrderIndex = 5;

        for l = 1:N_N
          hold on
          plot(log10_ks_combined{l}, omega_mu_k_combined{l}, 'LineWidth', 2)
          hold off
        end

        ax.ColorOrderIndex = 5;

        for l = 1:N_N
          hold on
          plot(log10_ks_combined{l}, omega_approx_k_combined{l}, '--', 'LineWidth', 2)
          hold off
        end

        xlabel('$$\log_{10}k$$', 'Interpreter', 'latex')
        ylabel('$$\frac{\mathrm{d}\Re\left(\omega\right)}{\mathrm{d}k}$$', 'Interpreter', 'latex')
        legend(legends_combined, 'Location', 'northeast', 'Interpreter', 'latex');
        title(title_combined_str, 'Interpreter', 'latex')
        axis tight
        ax.FontSize = 12;

        if save_flag == 1
          filename = ['review_ca_k', '_', save_str];
          saveas(gcf, [graphics_dir, filename, '.eps'], 'epsc')
        end

      end

    end

  end

end
