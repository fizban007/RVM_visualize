import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def cone_gaussian_intensity(r, theta, alpha, beta, phi, sigma, I0=1.0):
    # Cone axis rotated about z
    ax = np.sin(beta) * np.cos(phi)
    ay = np.sin(beta) * np.sin(phi)
    az = np.cos(beta)
    
    # Observation point (on x-z plane by construction)
    rx = r * np.sin(theta)
    ry = 0.0
    rz = r * np.cos(theta)
    term1 = np.abs((rx*rx + ry*ry + rz*rz) - (rx * ax + ry * ay + rz * az)**2)
    term2 = 0.2* 2 * np.tan(alpha+1e-10) **2
    I = I0 * np.exp(-term1/term2)
    
    return I


# Example parameters
r = 1.0
theta = np.pi / 3   # observation angle
alpha = np.pi / 6   # cone opening angle
# beta = np.pi / 4    # dipole angle
beta_range = np.linspace(0, np.pi/2, 10)
sigma = 0.1
I0 = 1.0

plt.figure(figsize=(6,4))
# Sweep phi from 0 to 2π
phi_vals = np.linspace(0, 4*np.pi, 400)
for beta in beta_range:
    intensity_vals = cone_gaussian_intensity(r, theta, alpha, beta, phi_vals, sigma, I0)
    plt.plot(phi_vals, intensity_vals, c=plt.cm.rainbow(beta/np.pi*2), label=f'Beta = {beta * 180/np.pi:.2f} deg')
# Plot
plt.xlabel(r"Rotation angle $\phi$ (rad)")
plt.ylabel("Intensity")
plt.title("Gaussian Cone Intensity vs Rotation")
plt.grid(True)
plt.legend()
plt.savefig("cone_intensity.png", dpi=300)
#plt.show()

