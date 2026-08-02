# spatial_SA_with_OT
Scripts and data of C. Perchtold's postdoc

- First load the groundwater sulphate concentration (file 01_spatial_points.Rdata) from aquifer 113AK in the Paris Basin, which has standardized coordinates.
- Then load the uncertainty matrix (file 01_uncertainty_matrix.Rdata), containing our uncertain factors: mesh, formula, covariance function, prior on the range and standard deviation. The generation of the meshes and formulas is not included.
- Step 1 to Step 4 explain our workflow applied to the case study (cross-validation with INLA, INLA, principal component analysis, sensitivity analysis with optimal transport) and should be used in this order. 

In case you do not want to run the simulations yourself, check out https://zenodo.org/records/21725111.
