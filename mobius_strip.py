import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def __init__(self, R, w, n):
        """
        Initialize MobiusStrip with radius R, width w, and resolution n.
        
        Args:
            R (float): Distance from center to strip midline
            w (float): Width of the strip
            n (int): Number of points in each dimension of the mesh
        """
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.mesh = self.compute_mesh()
        self.surface_area = self.compute_surface_area()
        self.edge_length = self.compute_edge_length()

    def parametric_equations(self, u, v):
        """
        Compute (x, y, z) coordinates for given u, v parameters.
        
        Args:
            u (ndarray): Parameter in [0, 2pi]
            v (ndarray): Parameter in [-w/2, w/2]
        
        Returns:
            tuple: (x, y, z) coordinates
        """
        x = (self.R + v * np.cos(u/2)) * np.cos(u)
        y = (self.R + v * np.cos(u/2)) * np.sin(u)
        z = v * np.sin(u/2)
        return x, y, z

    def compute_mesh(self):
        """
        Generate 3D mesh of points on the Mobius strip surface.
        
        Returns:
            tuple: (x, y, z) mesh arrays
        """
        return self.parametric_equations(self.U, self.V)

    def compute_surface_area(self):
        """
        Approximate surface area using numerical integration.
        Uses the cross product of partial derivatives to compute the area element.
        
        Returns:
            float: Approximated surface area
        """
        # Partial derivatives with respect to u and v
        def partial_derivatives(u, v):
            # du derivatives
            du_x = (-v * np.sin(u/2) * np.cos(u) / 2 
                   - (self.R + v * np.cos(u/2)) * np.sin(u))
            du_y = (-v * np.sin(u/2) * np.sin(u) / 2 
                   + (self.R + v * np.cos(u/2)) * np.cos(u))
            du_z = v * np.cos(u/2) / 2
            
            # dv derivatives
            dv_x = np.cos(u/2) * np.cos(u)
            dv_y = np.cos(u/2) * np.sin(u)
            dv_z = np.sin(u/2)
            
            return (du_x, du_y, du_z), (dv_x, dv_y, dv_z)
        
        # Compute cross product magnitude for area element
        du, dv = partial_derivatives(self.U, self.V)
        cross_product = np.array([
            du[1] * dv[2] - du[2] * dv[1],
            du[2] * dv[0] - du[0] * dv[2],
            du[0] * dv[1] - du[1] * dv[0]
        ])
        area_element = np.sqrt(np.sum(cross_product**2, axis=0))
        
        # Numerical integration using trapezoidal rule
        du = self.u[1] - self.u[0]
        dv = self.v[1] - self.v[0]
        return np.sum(area_element) * du * dv

    def compute_edge_length(self):
        """
        Compute the length of the Mobius strip's edge numerically.
        The edge consists of two boundaries at v = ±w/2.
        
        Returns:
            float: Total edge length
        """
        # Edge at v = w/2 and v = -w/2
        edge_length = 0
        for v in [self.w/2, -self.w/2]:
            x, y, z = self.parametric_equations(self.u, v)
            # Compute arc length by summing distances between consecutive points
            dx = np.diff(x)
            dy = np.diff(y)
            dz = np.diff(z)
            segment_lengths = np.sqrt(dx**2 + dy**2 + dz**2)
            edge_length += np.sum(segment_lengths)
        return edge_length

    def plot(self):
        """
        Visualize the Mobius strip in 3D using matplotlib.
        """
        x, y, z = self.mesh
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z, cmap='viridis', alpha=0.8)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Möbius Strip')
        plt.savefig('mobius_strip.png')
        plt.close()

def main():
    # Example usage
    R, w, n = 5.0, 2.0, 100
    mobius = MobiusStrip(R, w, n)
    print(f"Surface Area: {mobius.surface_area:.2f}")
    print(f"Edge Length: {mobius.edge_length:.2f}")
    mobius.plot()

if __name__ == "__main__":
    main()