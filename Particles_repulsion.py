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


def load_image(path):
    return Image.open(path).convert('L')  # Convert image to grayscale

@njit
def initialize_particles(map, num_particles):
    width, height = map.shape
    negative_coords = [(random.randint(0, width - 1), random.randint(0, height - 1)) for _ in range(num_particles)]
    to_be_removed = []
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            dist = np.sqrt((negative_coords[i][0] - negative_coords[j][0])**2 + (negative_coords[i][1] - negative_coords[j][1])**2)
            if dist < 1:
                to_be_removed.append(j)
    negative_coords = [p for i, p in enumerate(negative_coords) if i not in to_be_removed]

    # Add positive particles to the image with the density proportional to the "negative" pixel intensity
    positive_coords = np.zeros((width*height,3))
    # Coule be coarse grained if every > 1
    every = 1
    # These parameters could be adjusted for better representation
    charge_factor = 1
    charge_threshold = 0.1
    ppi=0
    for x in range(width//every):
        for y in range(height//every):
            intensity = map[x*every, y*every]
            charge = (1-intensity / 255) * charge_factor
            if charge > charge_threshold:
                # Particles could be distributed randomly within the square (pixel x every) x (pixel x every)
                # positive_coords[ppi] = np.array([x*every + np.random.randint(0, every), y*every + np.random.randint(0, every), charge])
                # or in a deterministic way
                positive_coords[ppi] = np.array([x*every, y*every, charge])
                ppi+=1
    positive_coords = positive_coords[:ppi]
    return negative_coords, positive_coords

@njit
def calculate_forces(Nparticles, Pparticles, cutoff):
    # Using numpy arrays for better compatibility and performance
    n = len(Nparticles)
    p = len(Pparticles)
    forces = np.zeros((n,2))

    # Repulsive forces between Nparticles
    for i in prange(n):
        for j in range(i + 1, n):
            dx = Nparticles[i][0] - Nparticles[j][0]
            dy = Nparticles[i][1] - Nparticles[j][1]
            distance = np.sqrt(dx**2 + dy**2) + 1e-5  # Prevent division by zero
            if distance < cutoff:
                force_magnitude = 1 / distance ** 2
                forces[i][0] += dx * force_magnitude
                forces[i][1] += dy * force_magnitude
                forces[j][0] -= dx * force_magnitude
                forces[j][1] -= dy * force_magnitude

    # Attractive or repulsive forces from Pparticles on Nparticles
    for i in prange(n):
        for j in range(p):
            dx, dy = Nparticles[i][0] - Pparticles[j][0], Nparticles[i][1] - Pparticles[j][1]
            distance = np.sqrt(dx**2 + dy**2) + 1e-5
            if distance < cutoff:
                force_magnitude = 1 / distance ** 2
                forces[i][0] -= dx * force_magnitude * Pparticles[j][2]
                forces[i][1] -= dy * force_magnitude * Pparticles[j][2]

    return forces

def move_particles(particles, forces,factor):
    new_positions = np.zeros((len(particles),2))
    for i,p in enumerate(particles):
        new_x = p[0] + forces[i][0]*factor
        new_y = p[1] + forces[i][1]*factor
        # new_x = max(0, min(new_x, img.width - 1))
        # new_y = max(0, min(new_y, img.height - 1))
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
    num_particles = 30000  
    cutoff = 20  # Cutoff distance for repulsive forces
    force_factor = 0.1
    total_frames = 25
    particle_size = 0.5
    DPI = 400
    Mesh = True
    # Save parameters in log file
    with open("log.txt", "w") as f:
        f.write(f"num_particles: {num_particles}\n")
        f.write(f"cutoff: {cutoff}\n")
        f.write(f"force_factor: {force_factor}\n")
        f.write(f"total_frames: {total_frames}\n")
        f.write(f"particle_size: {particle_size}\n")
        f.write(f"DPI: {DPI}\n")
        f.write(f"Mesh: {Mesh}\n")


    # Load image and show it
    img_path = sys.argv[1] 
    img = load_image(img_path)
    img.show()

    # Initialize particles
    map = np.array(img).transpose()
    Nparticles, Pparticles = initialize_particles(map, num_particles)

    fname = img_path.split('.')[0]
    def update(frame):
        nonlocal Nparticles
        forces = calculate_forces(Nparticles, Pparticles, cutoff)
        Nparticles = move_particles(Nparticles, forces, force_factor * Heaviside(frame))

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
            print(f"Frame {frame}/{total_frames} saved")
            return triplot
        else:            
            scat.set_offsets(Nparticles)        
            fig.canvas.draw()
            fig.savefig(f"{fname}_particle_frame_{frame+1:02d}.png", bbox_inches='tight', pad_inches=0, dpi=DPI)
            print(f"Frame {frame}/{total_frames} saved")
            return scat,

    # Initialize plot
    fig, ax = plt.subplots()
    if Mesh:
        # Construct Delaunay triangulation
        NPparticles = np.array([[p[0], p[1]] for p in Pparticles])
        tri = Delaunay(NPparticles)
        ax.triplot(NPparticles[:,0], NPparticles[:,1], tri.simplices, 'k-', linewidth=0.1)
    else:
        scat = ax.scatter(Nparticles[:][0], Nparticles[:][1], s=particle_size, c='black', linewidths=0, zorder=2)
    ax.set_xlim(0, img.width)
    ax.set_ylim(0, img.height)
    # set aspect ratio to respect the image size
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')
    # fig.savefig(f"{fname}_frame_00.png", bbox_inches='tight', pad_inches=0, dpi=200)

    # Create animation
    ani = FuncAnimation(fig, update, frames=total_frames, interval=200, blit=True)    
    if Mesh:
        ani.save(img_path.split('.')[0] + '_mesh_animation.mp4', writer='ffmpeg', dpi=200)
    else:
        ani.save(img_path.split('.')[0] + '_particle_animation.mp4', writer='ffmpeg', dpi=200)


if __name__ == "__main__":
    main()
