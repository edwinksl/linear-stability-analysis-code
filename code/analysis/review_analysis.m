function review_analysis(save_flag)

  close all

  % Define custom colormap
  addpath('DrosteEffect-BrewerMap-9e40245')
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Greens');
  map_3 = brewermap(9, 'Reds');
  map = [map_1([5,8],:); map_2([5,8],:); map_3(7,:)];

  bc_index = 1;
  base_dir = 'datasets/analysis/base/';
  lsa_dir = 'datasets/analysis/lsa/';
  graphics_dir = '../../review-graphics/';

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  rho_ss = [-0.05, 0, 0.05];
  Das = [0.1, 1, 10];
  J_applieds = [0.5, 1, 1.5];
  N = 1001;

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;
  omega_max_approx = max_approx_data.omega_max_approx;

  c_max_num_data = load([lsa_full_dir, 'c_max_num', '.mat']);
  k_c_num = c_max_num_data.k_c_num;
  k_max_num = c_max_num_data.k_max_num;
  omega_max_num = c_max_num_data.omega_max_num;

  % Collapse plots for numerical k_max, omega_max and k_c across rho_s values

  % Define custom colormaps
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Reds');
  map_k = map_1([4,6,8],:);
  map_omega = map_2([4,6,8],:);

  % line_styles = {'-', '--', ':'};

  for j = 1:N_Da

    Da = Das(j);

    for k = 1:N_J_applied

      J_applied = J_applieds(k);

      k_c_num_jk = k_c_num(:,j,k);
      k_max_num_jk = k_max_num(:,j,k);
      omega_max_num_jk = omega_max_num(:,j,k);

      % Plot V vs. t
      if save_flag == 1
        f_V = figure('Visible', 'off');
      else
        f_V = figure;
      end

      V_legend = cell(N_rho_s, 1);

      title_str = ['$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$J_\textnormal{a}$$ = ', num2str(J_applied, '%.3g')];  % rename J_applied to J_a in title

      filename_str = [num2str(Da), '_', num2str(J_applied), '_', num2str(N)];

      for i = 1:N_rho_s

        rho_s = rho_ss(i);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        filename = ['lsa', '_', str];
        lsa_data = load([lsa_full_dir, filename, '.mat']);

        base = lsa_data.base;
        t_indices = lsa_data.t_indices;

        t = base.t;
        t_s = base.t_s;
        rho_s = base.rho_s;
        Da = base.Da;
        J_applied = base.J_applied;

        c = base.c;
        phi = base.phi;
        c_1 = c(end, :);
        phi_1 = phi(end, :);
        eta_1 = base.eta_1;
        n = base.n;
        xi_p = base.xi_p;

        V = get_V(c_1, phi_1, eta_1, n, xi_p, rho_s);
        V_abs = abs(V);

        t_rs = t/t_s;  % rs = rescaled
        t_rs_num = t_rs(t_indices);

        V_legend{i} = ['$$\rho_\textnormal{s}$$ = ', num2str(rho_s)];

        % Plot V vs. t
        set(0, 'CurrentFigure', f_V)
        hold on
        plot(t_rs_num, V_abs, 'Color', map_omega(i,:), 'LineStyle', '-', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
        hold off

      end

      set(0, 'CurrentFigure', f_V)
      ax_V = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$\left|V\right|$$', 'Interpreter', 'latex')
      legend(V_legend, 'Location', 'southeast', 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      axis tight
      ax_V.FontSize = 14;

      if save_flag == 1
        V_filename = ['V', '_', filename_str];
        saveas(gcf, [graphics_dir, V_filename, '.eps'], 'epsc')
      end

    end

  end

end
