import numpy as np
import matplotlib.pyplot as plt

def cone_gaussian_intensity(r, xi, zeta, phi, sigma, I0=1.0):
    # Cone axis rotated about z
    ax = np.sin(xi) * np.cos(phi)
    ay = np.sin(xi) * np.sin(phi)
    az = np.cos(xi)

    # Observation point (on x-z plane by construction)
    rx = r * np.sin(zeta)
    ry = 0.0
    rz = r * np.cos(zeta)
    term1 = np.abs((rx*rx + ry*ry + rz*rz) - (rx * ax + ry * ay + rz * az)**2)
    term2 = 0.2* 2 * np.tan(sigma+1e-10) **2
    I = I0 * np.exp(-term1/term2)
    
    return I

def Polarization_Angle(zeta, phi, eta, eps, delta, alpha, beta):
    numer = (1 + eta - eps * np.cos(delta) * np.cos(zeta)) * np.sin(alpha) * np.sin(beta + phi) \
          + eps * np.sin(delta) * (np.cos(alpha) * np.cos(zeta) - np.sin(alpha) * np.sin(beta) * np.sin(zeta)) 
    denom = (1 + eta) * (np.cos(alpha) * np.sin(zeta) - np.sin(alpha) * np.cos(zeta) * np.cos(beta + phi)) \
          + eps * (np.sin(alpha) * np.cos(delta) * np.cos(beta + phi) - np.cos(alpha) * np.sin(delta) * np.cos(phi))
    return np.arctan2(numer,denom)


# Example parameters
r = 1.0
xi = np.pi / 3   # dipole angle
zeta_range = np.linspace(0, np.pi/2, 10) # angle between observer and spin axis
sigma = np.pi / 8 # cone opening angle
I0 = 1.0

alpha = xi	   #α   Magnetic moment obliquity	degrees
beta = np.deg2rad(40.0)	   #β   Magnetic moment phase at t=0	degrees
delta = np.deg2rad(0.0)	   #δ   Colatitude of translation vector	degrees
epsilon = 0.5  #ϵ   Normalized dipole translation distance (ϵ=d/R)	dimensionless
eta = 0.0	   #η   Normalized emission height (η=h/R)


phi_vals = np.linspace(0, 4*np.pi, 400)

# Create 1 column, 2 rows: intensity on top, PA on bottom
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 8))

# colormap for zeta
cmap = plt.cm.rainbow
norm = plt.Normalize(vmin=zeta_range.min(), vmax=zeta_range.max())
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

for zeta in zeta_range:
    intensity_vals = cone_gaussian_intensity(r, xi, zeta, phi_vals, sigma, I0)
    PA = Polarization_Angle(zeta, phi_vals, eta, epsilon, delta, alpha, beta)
    PA_deg = np.degrees(PA)
    color = cmap(norm(zeta))

    # plot intensity (top)
    ax1.plot(phi_vals, intensity_vals, color=color, alpha=0.9)

    # plot PA (bottom) with color where intensity>threshold, grey elsewhere
    # Use a small threshold since Gaussian intensity is never exactly 0
    threshold = 0.01 * I0
    
    # Create arrays for plotting with color masking
    PA_colored = PA_deg.copy()
    PA_grey = PA_deg.copy()
    
    # Mask: where intensity > threshold, use color; else use grey
    mask_color = intensity_vals > threshold
    mask_grey = ~mask_color
    
    # Set values to NaN where we don't want to plot them
    PA_colored[mask_grey] = np.nan
    PA_grey[mask_color] = np.nan
    
    # Plot colored segments (where intensity > threshold)
    ax2.plot(phi_vals, PA_colored, color=color, linestyle='-', alpha=0.9, linewidth=1.5)
    # Plot grey segments (where intensity <= threshold)
    ax2.plot(phi_vals, PA_grey, color='grey', linestyle='-', alpha=0.5, linewidth=1.0)

# labels and title
ax2.set_xlabel(r"Rotation angle $\phi$ (rad)")
ax1.set_xlim(phi_vals.min(), phi_vals.max())
ax1.set_ylabel("Intensity")
ax2.set_ylabel("PA (deg)")
ax2.set_ylim(-180, 180)
ax1.set_title(rf"Gaussian Cone Intensity and PA vs Rotation ($\xi = {xi * 180/np.pi:.1f}^\circ$)")

# colorbar showing zeta (in degrees) on the right
cbar = fig.colorbar(sm, ax=[ax1, ax2], pad=0.06, fraction=0.03)
ticks = np.linspace(zeta_range.min(), zeta_range.max(), 5)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{t*180/np.pi:.0f}°" for t in ticks])
cbar.set_label(r'$\zeta$ (deg)')

# legend explaining line meaning
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color='k', lw=2, label='Intensity (top)'),
                   Line2D([0], [0], color='k', lw=2, linestyle='-', label='PA (color where I>0)'),
                   Line2D([0], [0], color='grey', lw=2, linestyle='-', label='PA (I=0)')]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=8)

ax1.grid(True)
ax2.grid(True)
plt.savefig("cone_intensity.png", dpi=300)
# plt.show()

