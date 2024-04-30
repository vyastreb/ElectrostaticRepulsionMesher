"""
Code that converts images to halftone images using the particle repulsion algorithm.

Code based on the particle repulsion algorithm
[1] Schmaltz, C., Gwosdek, P., Bruhn, A. and Weickert, J., 2010, December. Electrostatic halftoning. In Computer Graphics Forum (Vol. 29, No. 8, pp. 2313-2327). Oxford, UK: Blackwell Publishing Ltd.
See also web-page: https://www.mia.uni-saarland.de/Research/Electrostatic_Halftoning/index.shtml

Author: Vladislav A. Yastrebov
Afffiliation: CNRS, MINES Paris - PSL
License: MIT
Date: April 26, 2024
Note: The code is provided as is and the author is not responsible for any damage caused by the code.
AI usage: GPT4 and copilot were used to co-construct the code.
"""
    
import numpy as np
from PIL import Image
import random
import sys
from numba import njit, prange
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.spatial import Delaunay
import datetime
    
def load_image(path):
    return Image.open(path).convert('L')  # Convert image to grayscale

@njit
def initialize_particles(map, Nx, Ny,factor_charge, zero_probability_offset, charge_threshold):
    width, height = map.shape
    negative_coords = np.zeros((Nx*Ny,2))
    positive_coords = np.zeros((Nx*Ny,3))
    dx = width/Nx
    dy = height/Ny
    nii = 0
    for i in range(Nx):
        x = dx*i 
        for j in range(Ny):
            y = dy*j
            if np.random.rand() > map[int(x), int(y)]/255 * (1 - zero_probability_offset): # + zero_probability_offset:
                # negative_coords[nii] = np.array([dx*i+np.random.rand()*dx, dy*j+np.random.rand()*dy])
                negative_coords[nii] = np.array([dx*(i+0.5), dy*(j+0.5)])
                nii += 1
    negative_coords = negative_coords[:nii]

    # Add positive particles to the image with the density proportional to the "negative" pixel intensity
    pii = 0
    total_charge = 0
    for i in range(Nx):
        for j in range(Ny):
            intensity = map[int(dx*i), int(dy*j)]
            charge = (1-intensity / 255)
            if charge > 1:
                print("Error: charge > 1")
            if charge > charge_threshold:
                # Particles could be distributed randomly within the square (pixel x every) x (pixel x every)
                positive_coords[pii] = np.array([dx*i + np.random.rand()*dx, dy*j + np.random.rand()*dy, charge])
                # or on a regular grid
                # positive_coords[pii] = np.array([dx*i, dy*j, charge])
                pii += 1
                total_charge += charge
    positive_coords = positive_coords[:pii]

    # To equilibrate the total charge of negative and positive particles
    factor = (nii / total_charge) / factor_charge
    positive_coords[:,2] *= factor

    return negative_coords, positive_coords

@njit
def calculate_forces(Nparticles, Pparticles, cutoff, eps_regulation):
    # Using numpy arrays for better compatibility and performance
    n = len(Nparticles)
    p = len(Pparticles)
    forces = np.zeros((n,2))

    # Repulsive forces between Nparticles
    for i in prange(n):
        for j in prange(i + 1, n):
            dx = Nparticles[i][0] - Nparticles[j][0]
            dy = Nparticles[i][1] - Nparticles[j][1]
            distance = np.sqrt(dx**2 + dy**2)  # Prevent division by zero
            if distance < cutoff:
                force_magnitude = 1 / (distance + eps_regulation) ** 2
                fx = dx * force_magnitude
                fy = dy * force_magnitude
                forces[i][0] += fx
                forces[i][1] += fy
                forces[j][0] -= fx
                forces[j][1] -= fy

    # Attractive or repulsive forces from Pparticles on Nparticles
    for i in prange(n):
        for j in prange(p):
            dx = Nparticles[i][0] - Pparticles[j][0]
            dy = Nparticles[i][1] - Pparticles[j][1]
            distance = np.sqrt(dx**2 + dy**2) 
            if distance < cutoff:
                force_magnitude = 1 / (distance + eps_regulation) ** 2
                forces[i][0] -= dx * force_magnitude * Pparticles[j][2]
                forces[i][1] -= dy * force_magnitude * Pparticles[j][2]

    return forces

def move_particles(img, particles, forces, factor, iter, max_iter):
    new_positions = np.zeros((len(particles),2))
    for i,p in enumerate(particles):
        normf = np.sqrt(forces[i][0]**2 + forces[i][1]**2)
        # Damping factor to prevent particles from moving back and forth around the same position, ad hoc solution. FIXME could be improved
        offset = max_iter - 10
        damping_factor = np.heaviside(iter - offset, 0) * (iter - offset) / 3. 

        new_x = p[0] + forces[i][0] * factor / (normf * (1 + damping_factor) )
        new_y = p[1] + forces[i][1] * factor / (normf * (1 + damping_factor) )

        # Prevent particles from moving outside the image
        new_x = max(1, min(new_x, img.width - 1))
        new_y = max(1, min(new_y, img.height - 1))

        new_positions[i] = np.array([new_x, new_y])
    return new_positions

def render_image(particles, img):
    fig,ax = plt.subplots()
    plt.scatter([p[0] for p in particles], [p[1] for p in particles], s=1)
    plt.xlim(0, img.width)
    plt.ylim(0, img.height)
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.show()

def Heaviside(x):
    if x > 0:
        return 1
    else:
        return 0

def main():
    # Parameters
    cutoff = 20  # Cutoff distance for electrostatic forces
    force_factor = 0.5
    total_frames = 20
    particle_size = 1.25
    DPI = 400
    eps_regulation = 0.1 # Regularization parameter to prevent division by zero (in pixel size)
    Mesh = True

    # Additional parameters, normally should be kept as is
    every = 1
    factor_charge           = 1.0
    zero_probability_offset = 0.02 # To keep white white, set to 0
    charge_threshold        = 0.0

    # Load image and show it
    img_path = sys.argv[1] 
    img = load_image(img_path)
    img.show()

    # Initialize particles
    map = np.array(img).transpose()
    Nparticles, Pparticles = initialize_particles(map, img.width//every, img.height//every, factor_charge, zero_probability_offset, charge_threshold)
    print("Number of negative particles:", Nparticles.shape[0])
    print("Number of positive particles:", Pparticles.shape[0])

    # Save parameters in log file
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    with open("log_"+timestamp+".txt", "w") as f:
        f.write(f"num_particles nagative: {Nparticles.shape[0]}\n")
        f.write(f"num_particles positive: {Pparticles.shape[0]}\n")
        f.write(f"eps_regulation: {eps_regulation}\n")
        f.write(f"cutoff: {cutoff}\n")
        f.write(f"force_factor: {force_factor}\n")
        f.write(f"total_frames: {total_frames}\n")
        f.write(f"particle_size: {particle_size}\n")
        f.write(f"zero_probability_offset: {zero_probability_offset}\n")
        f.write(f"DPI: {DPI}\n")
        f.write(f"Mesh: {Mesh}\n")
        f.write("--------------------\n")
        f.write(f"every (should be 1): {every}\n")
        f.write(f"charge_threshold (should be 0 to keep white white): {charge_threshold}\n")
        f.write(f"factor_charge (should be 1 to keep charge neutral box): {factor_charge}\n")


    # Nparticles, Pparticles = initialize_particles(map, num_particles)

    fname = img_path.split('.')[0]
    def update(frame):
        nonlocal Nparticles
        forces = calculate_forces(Nparticles, Pparticles, cutoff, eps_regulation)
        Nparticles = move_particles(img, Nparticles, forces, force_factor * Heaviside(frame), frame, total_frames)

        if Mesh:
            # Clear the plot and redraw the triangulation
            ax.clear()
            ax.set_xlim(0, img.width)
            ax.set_ylim(0, img.height)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.axis('off')
            # Recompute the Delaunay triangulation
            NPparticles = np.array([[p[0], p[1]] for p in Nparticles])
            tri = Delaunay(NPparticles)            
            triplot = ax.triplot(NPparticles[:,0], NPparticles[:,1], tri.simplices, 'k-', linewidth=0.1, zorder=3)            
            fig.canvas.draw()
            fig.savefig(f"{fname}_mesh_frame_{frame+1:02d}.png", bbox_inches='tight', pad_inches=0, dpi=DPI)
            print(f"Frame {frame+1}/{total_frames}")
            return triplot
        else:         
            # To show non-normalized forces uncomment these lines and comment `scat.set_offsets(Nparticles)`
            # ax.clear()   
            # ax.set_xlim(0, img.width)
            # ax.set_ylim(0, img.height)
            # ax.set_aspect('equal')
            # ax.invert_yaxis()
            # ax.axis('off')
            # ax.scatter(Nparticles[:,0], Nparticles[:,1], s=particle_size, c='black', linewidths=0, zorder=2)  
            # ax.quiver(Nparticles[:,0], Nparticles[:,1], forces[:,0], -forces[:,1], color='red', scale=10)
            scat.set_offsets(Nparticles)    
            fig.canvas.draw()            
            fig.savefig(f"{fname}_particle_frame_{frame+1:02d}.png", bbox_inches='tight', pad_inches=0, dpi=DPI)
            print(f"Frame {frame+1}/{total_frames}")
            return scat,

    # Initialize plot
    fig, ax = plt.subplots()
    # For transparent background uncomment these lines
    # fig.patch.set_alpha(0)
    # ax.patch.set_alpha(0)

    if Mesh:
        # Construct Delaunay triangulation
        NPparticles = np.array([[p[0], p[1]] for p in Pparticles])
        tri = Delaunay(NPparticles)
        ax.triplot(NPparticles[:,0], NPparticles[:,1], tri.simplices, 'k-', linewidth=0.1)
    else:
        scat = ax.scatter(Nparticles[:][0], Nparticles[:][1], s=particle_size, c='black', linewidths=0, zorder=2)
    ax.set_xlim(0, img.width)
    ax.set_ylim(0, img.height)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')

    # Create animation
    ani = FuncAnimation(fig, update, frames=total_frames, interval=200, blit=True)    
    if Mesh:
        ani.save(img_path.split('.')[0] + '_mesh_animation.mp4', writer='ffmpeg', dpi=200)
    else:
        ani.save(img_path.split('.')[0] + '_particle_animation.mp4', writer='ffmpeg', dpi=200)


if __name__ == "__main__":
    main()
