# Axial normal stress in an expanding transverse coordinate frame

Registered 2026-09-06 after mobile-source matching, before the new tests.
This is a coefficient-free straight-axis coordinate/product-identity module,
not a full curved Reynolds-stress or pressure closure. Do not modify the
frozen finite-TKE operator or relabel its omitted-normal approximation.

Use physical y_i=L_i(s)*xi_i, A proportional to L_y*L_n, beta_i=L_i'/L_i.
Write N=rho*Rss and T_i=rho*Rsi. The transformed fluxes are

    FM_i = A*rho*(v_i-u*L_i'*xi_i)/L_i
    FP_i = u*FM_i + A*T_i/L_i - A*N*beta_i*xi_i.

Thus FP-u*FM is not A*T/L once the axial normal stress is nonzero and the
width changes. Recover physical shear with the geometric term retained:

    T_i = L_i/A * (FP_i-u*FM_i + A*N*beta_i*xi_i).

Axial mean-kinetic-plus-normal-work flux is A*(rho*u^3/2+u*N).
Its relative transverse flux is u*FP-u^2*FM/2. The consistent axial-row
production is

    P = -(FP-u*FM).grad_xi(u) - A*N*(partial_s u at fixed xi)
      = -A*(T_i*partial_y_i u + N*partial_s u at fixed physical y).

The two coordinate forms must agree. For arbitrary smooth fields (without
assuming their mass/momentum residuals vanish), the product identity is

    div(KE_flux) + P = u*div(momentum_flux) - u^2*div(mass_flux)/2.

Implement forward transformation, inverse physical shear/covariance and the
production term with strict finite/shape/positive-area-length-density guards.
Require explicit straight-axis geometry, rejecting curvature rather than
silently dropping its metric terms. Do not clip negative production.

Test independent physical polynomial fields and exponentially changing widths.
Evaluate actual flux divergences by central differences in s,xi_y,xi_n at
h=1e-4 and5e-5. Require scaled product error<=1e-7 at both steps and the
two production expressions agree to1e-12. Test inverse shear recovery, zero
normal-stress and zero-width-rate limits, and physical-coordinate velocity
gradient recovery. Include nonzero mass/momentum residuals in the identity.

The downstream operator must still be reassembled with N derivatives,
pressure, transverse momentum/stress work, curvature, ambient stress data
and a declared turbulence closure. These helpers return no trajectory rates,
physical acceptance or observation score. No synthetic stress becomes an inlet.
