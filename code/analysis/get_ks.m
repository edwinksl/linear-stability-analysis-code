function ks = get_ks(n_min, n_max, delta_n, L_y, L_z, k_cutoff)

  n_cutoff = round(get_n_approx(k_cutoff, L_y, L_z));

  n_range = [n_min:1:n_cutoff, n_cutoff+1:delta_n:n_max];  % use delta_n = 1 first
  n_y = n_range;
  n_z = n_range;
  ks = sqrt((n_y*pi/L_y).^2 + (n_z*pi/L_z).^2);

end
