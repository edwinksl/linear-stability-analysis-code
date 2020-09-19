function copper_analysis(save_flag)

  close all

  copper_dir = 'datasets/copper/';
  graphics_dir = '../../graphics/';

  data = load([copper_dir, 'copper.mat']);
  data = data.data;

  h_c = data.h_c;
  h_c_lb = data.h_c_lb;
  h_c_ub = data.h_c_ub;
  I_dims = data.I_dims;
  lambda_c_dims = data.lambda_c_dims;
  lambda_max_dims = data.lambda_max_dims;
  omega_max_dims = data.omega_max_dims;

  h_c = h_c*1e6;  % convert m to um
  h_c_lb = h_c_lb*1e6;  % convert m to um
  h_c_ub = h_c_ub*1e6;  % convert m to um
  I_dims = I_dims*1e3;  % convert A to mA
  lambda_c_dims = lambda_c_dims*1e6;  % convert m to um
  lambda_max_dims = lambda_max_dims*1e6;  % convert m to um

  I_dim_min = min(I_dims);
  I_dim_max = max(I_dims);

  lambda_c_min = 0;
  lambda_c_max = max(max(lambda_c_dims));
  lambda_c_max_rounded = 1.3;

  lambda_max_min = 0;
  lambda_max_max = max(max(lambda_max_dims));
  lambda_max_max_rounded = 6.4;

  omega_max_min = 0;
  omega_max_dims_finite = omega_max_dims(isfinite(omega_max_dims));
  omega_max_max = max(max(omega_max_dims_finite));
  omega_max_max_rounded = 1e-2;

  if save_flag == 1
    figure('Visible', 'off')
  else
    figure
  end

  patch([I_dim_min, I_dim_max, I_dim_max, I_dim_min], [h_c, h_c, lambda_c_max_rounded, lambda_c_max_rounded], 'blue', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')
  patch([I_dim_min, I_dim_max, I_dim_max, I_dim_min], [h_c, h_c, lambda_c_min, lambda_c_min], 'red', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')

  ax = gca;
  hold on
  plot(I_dims, lambda_c_dims, '-o', 'LineWidth', 2)
  plot(I_dims, h_c*ones(size(I_dims)), 'k--', 'LineWidth', 2)
  lb = plot(I_dims, h_c_lb*ones(size(I_dims)), 'k--', 'LineWidth', 2);
  ub = plot(I_dims, h_c_ub*ones(size(I_dims)), 'k--', 'LineWidth', 2);
  lb.Color(4) = 0.25;
  ub.Color(4) = 0.25;
  hold off
  xlabel('$$I_\textnormal{a}\,/\,\textnormal{mA}$$', 'Interpreter', 'latex')
  ylabel('$$\lambda_\textnormal{c}\,/\,\mu\textnormal{m}$$', 'Interpreter', 'latex')
  legend({'CN(-)', 'CN(+)', '$$h_\textnormal{c}$$ with LB and UB'}, 'Interpreter', 'latex')
  axis([I_dim_min, I_dim_max, lambda_c_min, lambda_c_max_rounded])
  ax.FontSize = 14;
  ax.TickDir = 'out';  % tick marks directed inward are not visible in EPS files for unknown reasons

  if save_flag == 1
    copper_lambda_c_filename = 'copper_lambda_c';
    saveas(gcf, [graphics_dir, copper_lambda_c_filename, '.eps'], 'epsc')
  end

  if save_flag == 1
    figure('Visible', 'off')
  else
    figure
  end

  ax = gca;
  plot(I_dims, lambda_max_dims, '-o', 'LineWidth', 2)
  xlabel('$$I_\textnormal{a}\,/\,\textnormal{mA}$$', 'Interpreter', 'latex')
  ylabel('$$\lambda_\textnormal{max}\,/\,\mu\textnormal{m}$$', 'Interpreter', 'latex')
  legend({'CN(-)', 'CN(+)'})
  axis([I_dim_min, I_dim_max, lambda_max_min, lambda_max_max_rounded])
  ax.FontSize = 14;

  if save_flag == 1
    copper_lambda_max_filename = 'copper_lambda_max';
    saveas(gcf, [graphics_dir, copper_lambda_max_filename, '.eps'], 'epsc')
  end

  if save_flag == 1
    figure('Visible', 'off')
  else
    figure
  end

  ax = gca;
  plot(I_dims, omega_max_dims, '-o', 'LineWidth', 2)
  xlabel('$$I_\textnormal{a}\,/\,\textnormal{mA}$$', 'Interpreter', 'latex')
  ylabel('$$\omega_\textnormal{max}\,/\,\textnormal{s}^{-1}$$', 'Interpreter', 'latex')
  legend({'CN(-)', 'CN(+)'}, 'Location', 'northwest')
  axis([I_dim_min, I_dim_max, omega_max_min, omega_max_max_rounded])
  ax.YAxis.Exponent = -2;
  ax.FontSize = 14;

  if save_flag == 1
    copper_omega_max_filename = 'copper_omega_max';
    saveas(gcf, [graphics_dir, copper_omega_max_filename, '.eps'], 'epsc')
  end

end
