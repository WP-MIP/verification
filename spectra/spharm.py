def do_transform_scalar(field, ntrunc=None, assume_latlon=None, flip_lat=False):
    import numpy as np
    import pyshtools as pysh
    """
    Compute total power (En) and power per mode (Pn = En/(2l+1))
    from a scalar field on a global Gaussian grid using GLQ and pyshtools.
    """
    field = np.asarray(field)
    ny, nx = field.shape

    # 1. Determine dimension ordering
    # -------------------------------
    # The input field might be given as (lat, lon) or (lon, lat), and we want
    # to work internally as (nlat, nlon) = (latitude index, longitude index).
    #
    # - If the user explicitly tells us the orientation via `assume_latlon`,
    #   we trust that.
    # - Otherwise, we infer it by assuming that the *smaller* dimension is
    #   latitude (common for global grids: Nlat < Nlon).
    if assume_latlon == 'latlon':
        # Field is already organized as (lat, lon)
        nlat, nlon = ny, nx
    elif assume_latlon == 'lonlat':
        # Field is (lon, lat) → transpose to (lat, lon)
        nlat, nlon = nx, ny
        field = field.T
    else:
        # Infer: smaller dimension is latitude
        if ny <= nx:
            # Already (lat, lon)
            nlat, nlon = ny, nx
        else:
            # Assume (lon, lat) and transpose
            nlat, nlon = nx, ny
            field = field.T

    # Optionally flip latitude (based on how Gaussian grid is defined in input file)
    # Many model output conventions store latitude from North to South; some tools
    # (like pyshtools GLQ) expect South to North. `flip_lat=True` allows handling
    # this convention difference.
    if flip_lat:
        field = field[::-1, :]

    # 2. Set lmax and enforce GLQ shape (nlat = lmax+1, nlon = 2*lmax+1)
    # -------------------------------------------------------------------
    # For Gaussian-Legendre Quadrature (GLQ) used by pyshtools:
    #   - The number of Gaussian latitudes is nlat = lmax + 1
    #   - The number of longitudes required is nlon = 2*lmax + 1
    #
    # From the given grid, we deduce lmax and check that nlon is compatible.
    lmax = nlat - 1
    required_nlon = 2 * lmax + 1

    # If more longitudes than required, truncate in longitude.
    # If fewer, the grid is not suitable for SHGLQ and we abort.
    if nlon > required_nlon:
        # Drop extra longitudes to match SHGLQ requirements
        field_glq = field[:, :required_nlon]
    elif nlon < required_nlon:
        # Not enough longitudes for accurate GLQ-based spherical transform
        raise ValueError(
            f"Grid not sufficient for GLQ: need nlon={required_nlon}, have {nlon}."
        )
    else:
        # Exactly the required number of longitudes
        field_glq = field

    # Ensure the array is contiguous in memory (required by some Fortran-based
    # routines for efficiency and correctness).
    field_glq = np.ascontiguousarray(field_glq)

    # Set effective truncation:
    # - If user did not specify ntrunc, use the maximum allowed by the grid (lmax).
    # - Otherwise, keep the smaller of (user ntrunc, lmax).
    if ntrunc is None:
        ntrunc = lmax
    ntrunc = min(ntrunc, lmax)

    # 3. GLQ nodes/weights and spherical-harmonic expansion
    # -----------------------------------------------------
    # SHGLQ(lmax) returns:
    #   x : 1D array of Gaussian nodes = sin(latitude)
    #   w : corresponding quadrature weights
    #
    # We then perform the scalar spherical harmonic transform using SHExpandGLQ.
    x, w = pysh.expand.SHGLQ(lmax)

    # cilm shape: (2, lmax+1, lmax+1)
    # --------------------------------
    # cilm[0, l, m] = C_{l,m}, real cosine coefficients
    # cilm[1, l, m] = S_{l,m}, real sine coefficients
    # For each degree l: m = 0..l are defined, m > l are unused (remain zero).
    cilm = pysh.expand.SHExpandGLQ(field_glq, w, x, lmax_calc=lmax)

    # 4. Compute spectrum: total power En and power per mode Pn
    # ---------------------------------------------------------
    # Total power at degree l:
    #   En[l] = sum_{m=0}^l (C_{l,m}^2 + S_{l,m}^2)
    #
    # The "per-mode" spectrum Pn is En divided by the number of modes at that
    # degree, which is (2l+1). This is useful for comparing energy density
    # per mode instead of total energy per degree.
    En = np.zeros(ntrunc + 1)

    for l in range(ntrunc + 1):
        # Truncate coefficients up to m = l for this degree
        C_l = cilm[0, l, :l+1]
        S_l = cilm[1, l, :l+1]
        # Sum of squared magnitudes of (C_{l,m}, S_{l,m}) over all m
        En[l] = np.sum(C_l**2 + S_l**2)

    # Degree array (total wavenumber index)
    degrees = np.arange(ntrunc + 1)

    return degrees, En



def get_spectral_diff_matrix(x, n):
    import pyshtools as pysh
    import numpy as np
    """
    Constructs the Gauss-Legendre differentiation matrix D.
    
    The matrix D maps function values at Gauss-Legendre nodes (roots of P_n)
    to the values of the derivative at those same nodes.
    
    Parameters
    ----------
    x : np.ndarray
        Array of n Gauss-Legendre nodes (roots of P_n).
    n : int
        The number of nodes (degree of the grid).
        
    Returns
    -------
    D : np.ndarray
        The (n x n) differentiation matrix.
    """
    
    # 1. Calculate P_{n-1}(x) at the nodes.
    # We need P'_n(x) for the formula. Using the recurrence relation:
    # (1 - x^2) * P'_n(x) = n * P_{n-1}(x) - n * x * P_n(x)
    # Since x are roots of P_n, P_n(x) = 0. Thus: P'_n(x) = n * P_{n-1}(x) / (1 - x^2).
    
    # pysh.legendre.PLegendre(L, z) returns [P_0(z), ..., P_L(z)].
    # We request order n-1 and take the last element [-1] to get P_{n-1}.
    pp = np.array([pysh.legendre.PLegendre(n - 1, val)[n - 1] for val in x])
    
    # 2. Compute P'_n(x) at the nodes using the recurrence derived above.
    # Note: 1 - x^2 is the denominator factor.
    p_prime = (n * pp) / (1 - x**2)
    
    # 3. Construct the Off-Diagonal elements (i != j).
    # Formula: D_ij = (P'_n(x_i) / P'_n(x_j)) * (1 / (x_i - x_j))
    
    # Create meshgrids to compute all pairs (x_i - x_j)
    xi, xj = np.meshgrid(x, x, indexing='ij')
    dx = xi - xj # matrix of coordinate distances (x_i-x_j)
    
    # Add identity to dx to avoid division by zero on diagonal (temporarily)
    # D_ij = col_vector / row_vector / diff_matrix
    # Here, col_vector/row_vector creates the matrix of ratios: P'_n(x_i)/P'_n(x_j)
    D = (p_prime[:, np.newaxis] / p_prime[np.newaxis, :]) / (dx + np.eye(n))
    
    # 4. Handle Diagonal elements (i == j).
    # Analytically, D_ii = x_i / (1 - x_i^2).
    # However, for numerical stability and to ensure the derivative of a constant
    # is exactly zero, we use the identity: sum_j(D_ij) = 0.
    # Therefore, D_ii = - sum_{j != i} D_ij.
    
    # First, zero out the diagonal which currently contains garbage from the division step
    np.fill_diagonal(D, 0)
    
    # Replace diagonal elements with the negative sum of the off-diagonal elements of the row
    D[np.diag_indices_from(D)] = -np.sum(D, axis=1)
    
    return D



def compute_vorticity_divergence_on_grid (u, v, nlat, required_nlon, a, phi, cosphi, w, x):
    import numpy as np
   
    # Compute vorticity and divergence on the grid
    # -----------------------------------------------
    # We first compute the horizontal vorticity (ζ) and divergence (δ) fields
    # from U and V on the Gaussian grid:
    #
    #   λ = longitude, φ = latitude
    #   ζ = (1 / (a cosφ)) * ( ∂v/∂λ - ∂(u cosφ)/∂φ )
    #   δ = (1 / (a cosφ)) * ( ∂u/∂λ + ∂(v cosφ)/∂φ )
    #
    # We compute:
    #   - Longitudinal derivatives ∂/∂λ via FFT in longitude
    #   - Latitudinal derivatives ∂/∂φ via centered finite differences.
    vort = np.zeros_like(u)
    div  = np.zeros_like(u)

    D=get_spectral_diff_matrix(x,nlat)

    # Set up wave numbers for FFT in longitude.
    # np.fft.fftfreq(n, d) returns frequencies in cycles per unit of `d`.
    # Here d = 1.0/required_nlon → grid step in index space; the factor 2πi
    # is not included, so we use `1j * k` to represent ∂/∂λ in spectral space.
    ik = 1j * np.fft.fftfreq(required_nlon, d=1.0 / required_nlon)

    # Loop over latitudes and compute longitudinal derivatives via FFT
    for j in range(nlat):
        # Take FFT along longitude, multiply by ik to get derivative in spectral
        # space, then inverse FFT to get derivative in physical space.
        du_dlam = np.fft.ifft(ik * np.fft.fft(u[j, :])).real
        dv_dlam = np.fft.ifft(ik * np.fft.fft(v[j, :])).real

        # First terms in ζ and δ formulas (1/(a cosφ) factor):
        vort[j, :] = dv_dlam / (a * cosphi[j])
        div[j, :]  = du_dlam / (a * cosphi[j])

    # ∂/∂φ vectorized for all longitudes:
    # shape: (nlat, nlon)
    f_u = u * cosphi[:, None]
    f_v = v * cosphi[:, None]

    d_ucosphi_dphi = D @ f_u   # (nlat,nlat) @ (nlat,nlon) -> (nlat,nlon)
    d_vcosphi_dphi = D @ f_v
    

    vort -= d_ucosphi_dphi / a
    div  += d_vcosphi_dphi / a

    return vort, div




def do_transform_ke(field_u, field_v, ntrunc=None, assume_latlon=None, flip_lat=False):
    import pyshtools as pysh
    import numpy as np
    """
    Compute isotropic kinetic energy spectrum from U, V winds
    on a Gaussian grid using spherical harmonics (pyshtools).

    Input:
      field_u, field_v : arrays (lon, lat)
      ntrunc           : spectral truncation (optional)

    Output:
      degrees : total wavenumber n
      En_tot  : total kinetic energy spectrum
      En_rot  : rotational (vortical) KE spectrum
      En_div  : divergent KE spectrum
    """
    u = np.asarray(field_u)
    ny, nx = u.shape
    
    # 1. Determine dimension ordering
    # -------------------------------
    # Similar logic as in the scalar case:
    #   - Make sure U and V are organized as (lat, lon) internally.
    #   - Either use the explicit `assume_latlon` flag or infer from dimensions.
    if assume_latlon == 'latlon':
        # Fields are already (lat, lon)
        nlat, nlon = ny, nx
    elif assume_latlon == 'lonlat':
        # Fields are (lon, lat) → transpose both U and V
        nlat, nlon = nx, ny
        u = np.asarray(field_u).T
        v = np.asarray(field_v).T
    else:
        # Infer ordering by assuming the smaller dimension is latitude.
        if ny <= nx:
            # (lat, lon)
            nlat, nlon = ny, nx
        else:
            # (lon, lat) → transpose both U and V
            nlat, nlon = nx, ny
            u = np.asarray(field_u).T
            v = np.asarray(field_v).T

    # Optionally flip latitude (based on how Gaussian grid is defined in input file)
    # If data is stored from North→South but GLQ expects South→North, we flip.
    if flip_lat:
        u = u[::-1, :]
        v = v[::-1, :]

    # 2. Set lmax and enforce GLQ shape (nlat = lmax+1, nlon = 2*lmax+1)
    # -------------------------------------------------------------------
    # As in the scalar case:
    #   lmax = nlat - 1
    #   required_nlon = 2*lmax + 1
    #
    # Adjust or reject based on longitude count.
    lmax = nlat - 1
    required_nlon = 2 * lmax + 1

    if nlon > required_nlon:
        # Truncate extra longitudes to match GLQ requirements
        u = u[:, :required_nlon]
        v = v[:, :required_nlon]
    elif nlon < required_nlon:
        # Not enough longitudes for GLQ-based transforms
        raise ValueError(
            f"Grid not sufficient for GLQ: need nlon={required_nlon}, have {nlon}"
        )

    # Ensure memory contiguity for Fortran-style routines
    u = np.ascontiguousarray(u)
    v = np.ascontiguousarray(v)

    # Determine spectral truncation
    if ntrunc is None:
        ntrunc = lmax
    ntrunc = min(ntrunc, lmax)


    #  3. GLQ nodes/weights and spherical-harmonic expansion setup
    # ------------------------------------------------------------
    # Get GLQ nodes (x) and weights (w) for degree up to lmax.
    x, w = pysh.expand.SHGLQ(lmax)
    # GLQ node x is sin(latitude). Invert to latitude φ (radians).
    phi = np.arcsin(x)          # latitude [rad]
    # Precompute cos(latitude), used repeatedly for metric factors.
    cosphi = np.cos(phi)

    # Earth's radius [m] used for metric terms and KE scaling
    a = 6.37122e6                 # Earth radius [m]


    # 4. Compute vorticity and divergence on the grid
    # -----------------------------------------------
    vort, div = compute_vorticity_divergence_on_grid (
            u, v, nlat, required_nlon, a, phi, cosphi, w, x)

    # 5. Spherical harmonic transforms of vorticity and divergence
    # ------------------------------------------------------------
    # Expand ζ and δ into spherical harmonics:
    #   vort_cilm, div_cilm have shape (2, lmax+1, lmax+1),
    #   where [:, l, m] store (C_{l,m}, S_{l,m}) coefficients.
    vort_cilm = pysh.expand.SHExpandGLQ(vort, w, x, lmax_calc=lmax)
    div_cilm  = pysh.expand.SHExpandGLQ(div,  w, x, lmax_calc=lmax)

    # 6. Invert Laplacian → streamfunction (ψ) and velocity potential (χ)
    # --------------------------------------------------------------------
    # On the sphere, vorticity and divergence are related to the streamfunction
    # ψ and velocity potential χ by the spherical Laplacian:
    #
    #   ∇²ψ = ζ
    #   ∇²χ = δ
    #
    # In spherical harmonics, the Laplacian acting on degree l gives:
    #   ∇² Y_{l,m} = -l(l+1)/a² * Y_{l,m}
    #
    # Thus, in spectral space,
    #   ψ_{l,m} = -a² / [l(l+1)] * ζ_{l,m}
    #   χ_{l,m} = -a² / [l(l+1)] * δ_{l,m}
    #
    # We skip l=0 (no meaning for streamfunction/potential at l=0).
    psi = np.zeros_like(vort_cilm)
    chi = np.zeros_like(div_cilm)

    for l in range(1, lmax + 1):
        # Factor corresponding to inversion of Laplacian for degree l
        factor = -a**2 / (l * (l + 1))
        # Apply factor to both C and S coefficients up to order m = l
        psi[:, l, :l+1] = factor * vort_cilm[:, l, :l+1]
        chi[:, l, :l+1] = factor * div_cilm[:, l, :l+1]

    # 7. KE spectra (total, rotational, divergent)
    # -------------------------------------------
    # Kinetic energy can be decomposed into:
    #   - Rotational (vortical) KE, associated with streamfunction ψ
    #   - Divergent KE, associated with velocity potential χ
    #
    # For each total wavenumber (degree) l:
    #   E_rot(l) = [l(l+1) / (2 a²)] * Σ_m |ψ_{l,m}|²
    #   E_div(l) = [l(l+1) / (2 a²)] * Σ_m |χ_{l,m}|²
    #   E_tot(l) = E_rot(l) + E_div(l)
    #
    # We start from l=1 since l=0 contains no KE in this decomposition.
    En_tot = np.zeros(ntrunc + 1)
    En_rot = np.zeros(ntrunc + 1)
    En_div = np.zeros(ntrunc + 1)

    for l in range(1, ntrunc + 1):
        # Common factor in KE spectrum formulas at degree l
        prefac = l * (l + 1) / (2.0 * a**2)

        # Sum of squared streamfunction and velocity potential coefficients
        # over m = 0..l for both cosine (index 0) and sine (index 1) parts.
        psi_sq = np.sum(psi[0, l, :l+1]**2 + psi[1, l, :l+1]**2)
        chi_sq = np.sum(chi[0, l, :l+1]**2 + chi[1, l, :l+1]**2)

        # Rotational and divergent KE at degree l
        En_rot[l] = prefac * psi_sq
        En_div[l] = prefac * chi_sq
        En_tot[l] = En_rot[l] + En_div[l]

    # Degree array (total wavenumber index)
    degrees = np.arange(ntrunc + 1)

    return degrees, En_tot, En_rot, En_div

