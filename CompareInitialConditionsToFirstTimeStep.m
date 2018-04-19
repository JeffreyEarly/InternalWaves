ic_file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s/input/SaveIC_EarlyIWmodel.nc';
t0_file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s/output/3D/CONCAT/XYZ_000000.nc';

u_ic = ncread(ic_file,'u');
v_ic = ncread(ic_file,'v');
rho_prime_ic = ncread(ic_file,'s1');
rho_bar_ic = ncread(ic_file,'rhobar');
rho_ic = rho_prime_ic + reshape(rho_bar_ic,1,1,[]);


u_t0 = ncread(t0_file,'u');
v_t0 = ncread(t0_file,'v');
rho_t0 = ncread(t0_file,'rho');
rho_prime_t0 = rho_t0 - reshape(rho_bar_ic,1,1,[]);

error = @(u,u_unit) max(max(max( abs(u-u_unit)/ max(max(max(abs(u_unit))))   )));
error(u_ic,u_t0)
error(v_ic,v_t0)
error(rho_ic,rho_t0)
error(rho_prime_ic,rho_prime_t0)