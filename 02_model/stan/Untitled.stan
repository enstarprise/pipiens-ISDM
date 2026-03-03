data {
  int n_land_covs;
  int n_grids_total;
  matrix[n_grids_total, n_land_covs] z_land;
}
parameters {
  vector[n_land_covs] beta_land;
}
model {
  beta_land ~ normal(0, 0.1);
  
  // Debug
  if (iteration() == 1) {
    print("Test: beta_land = ", beta_land);
    print("z_land[1,] = ", z_land[1,]);
    print("dot_product = ", dot_product(beta_land, z_land[1,]));
  }
}