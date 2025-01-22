"""
Arc path planner
"""

import numpy as np


def arc(x0, u, N, dt):
    """Generates an arc trajectory given an initial state and control input"""

    def dx(v, w, t, eps=1e-15):
        return v * np.sin((w + eps) * t) / (w + eps)

    def dy(v, w, t, eps=1e-15):
        return -v * (1 - np.cos((w + eps) * t)) / (w + eps)

    traj = np.zeros((N, 3))
    traj[:, 0] = [x0[0] + dx(u[0], u[1], i * dt) for i in range(N)]
    traj[:, 1] = [x0[1] + dy(u[0], u[1], i * dt) for i in range(N)]
    traj[:, 2] = x0[2] + u[1] * np.arange(N) * dt

    return traj


class ArcPathPlanner:
    """Local path planner that selects the lowest cost path from a set of candidate arcs

    Takes a costmap and goal point as input, both assumed to be in local coordinates.
    Currently assumes constant linear speed, and thus boils down to selecting and angular velocity for the arc.

    """

    def __init__(
        self,
        dem,
        speed: float = 3.0,  # m/s
        max_omega: float = 0.3,  # rad/s
        arc_duration: float = 5.0,  # seconds
        dt: float = 0.1,  # seconds, used to sample points along the arc for cost evaluation
        num_arcs: int = 15,  # number of candidate arcs to consider
        costmap_shape: tuple = (31, 31),  # dimensions of the costmap
        costmap_resolution: float = 1.0,  # meters per pixel
    ):
        # Arc parameters
        self.dem = dem
        self.arc_duration = arc_duration  # seconds
        self.num_arc_points = int(self.arc_duration / dt)
        self.speed = speed
        self.omegas = np.linspace(-max_omega, max_omega, num_arcs)
        self.candidate_arcs = [
            arc(np.zeros(3), [self.speed, w], self.num_arc_points, dt)
            for w in self.omegas
        ]

        # Costmap parameters
        self.costmap_shape = costmap_shape
        self.costmap_resolution = costmap_resolution
        self.costmap_center = np.array(costmap_shape) / 2
        self.costmap = np.zeros(costmap_shape)

        self.goal = np.zeros(2)

    def update_goal(self, goal: np.ndarray):
        self.goal = goal

    def update_costmap(self, costmap: np.ndarray):
        self.costmap = costmap
        # print(self.costmap)

    # def update_costmap_from_dem(self, dem: np.ndarray):
    #     # TODO: Compute costmap from DEM based on slope
    #     resolution = 1.0
    #     gradients_x, gradients_y = np.gradient(dem, resolution)
    #     slope = np.sqrt(gradients_x**2 + gradients_y**2)
    #     roughness = np.abs(gradients_x - gradients_y)

    #         # Combine slope and roughness into a single cost
    #     costmap = 10 * slope + 5 * roughness  # Weight slope and roughness
    #     costmap = np.clip(costmap, 0, 255)    # Normalize costs for visualization

    #     return costmap

    #     self.costmap = dem

    def costmap_val(self, x: float, y: float):
        """Returns the value of the costmap at the given position, rounding to the nearest index."""
        cx, cy = self.costmap_center
        if isinstance(x, np.ndarray):
            # Use np.round to get the closest indices
            x_idx = np.round(-x + cx).astype(int)
            y_idx = np.round(-y + cy).astype(int)
            return self.costmap[x_idx, y_idx]
        else:  # Scalar values
            x_idx = int(np.round(-x + cx))
            y_idx = int(np.round(-y + cy))
            return self.costmap[x_idx, y_idx]

    def find_costmap(self, center, heading_angle):
        dem = self.dem
        resolution = self.costmap_resolution
        cx, cy = center
        size = self.costmap_shape[0]
        half_size = (size - 1) / 2

        print(center, size, heading_angle, resolution)
        # Define the box corners relative to the center
        local_corners = np.array(
            [
                [-half_size, -half_size],  # Bottom left
                [half_size, -half_size],  # Bottom right
                [half_size, half_size],  # Top right
                [-half_size, half_size],  # Top left
            ]
        )

        # Rotation matrix to align the box with the heading angle
        rotation_matrix = np.array(
            [
                [np.cos(heading_angle), -np.sin(heading_angle)],
                [np.sin(heading_angle), np.cos(heading_angle)],
            ]
        )

        # Rotate and translate the corners
        rotated_corners = (rotation_matrix @ local_corners.T).T + np.array([cx, cy])

        # Create a grid in the local frame
        x_coords = np.arange(-half_size, half_size + resolution, resolution)
        y_coords = np.arange(-half_size, half_size + resolution, resolution)
        grid_x, grid_y = np.meshgrid(x_coords, y_coords)
        local_grid_points = np.vstack([grid_x.ravel(), grid_y.ravel()]).T
        rotated_grid_points = (rotation_matrix @ local_grid_points.T).T + np.array(
            [cx, cy]
        )

        # Compute DEM indices by rounding the rotated grid points
        grid_indices = np.round(rotated_grid_points).astype(int)

        gradients_x, gradients_y = np.gradient(dem)
        gradient_magnitudes = np.sqrt(gradients_x**2 + gradients_y**2)

        gradient_costmap = np.full((size, size), np.nan)
        dem_costmap = np.full((size, size), np.nan)

        # Extract a 2D costmap from the DEM for the rotated grid
        dem_costmap = np.full((size, size), np.nan)
        for idx, (y, x) in enumerate(grid_indices):
            row = idx // size
            col = idx % size
            if 0 <= y < dem.shape[0] and 0 <= x < dem.shape[1]:
                gradient_costmap[row, col] = gradient_magnitudes[y, x]
                dem_costmap[row, col] = dem[y, x]
            # else:
            #     print(f"Out-of-bounds grid index: ({y}, {x})")

        return (
            rotated_corners,
            grid_indices,
            dem_costmap,
            gradient_costmap,
            rotation_matrix,
        )

    def plan(self):
        """Selects the best arc based on cost evaluation"""
        costs = np.zeros(len(self.omegas))
        costmap_costs = np.zeros(len(self.omegas))
        goal_costs = np.zeros(len(self.omegas))
        for i, w in enumerate(self.omegas):
            arc = self.candidate_arcs[i]
            # Steering cost
            costs[i] += np.abs(w)
            # Costmap cost
            cmap_cost = (
                np.sum(self.costmap_val(arc[:, 0], arc[:, 1])) / self.num_arc_points
            )
            costs[i] += cmap_cost
            costmap_costs[i] = cmap_cost
            # Goal cost
            goal_cost = np.linalg.norm(self.goal - arc[-1, :2])
            costs[i] += goal_cost
            goal_costs[i] = goal_cost

        idx = np.argmin(costs)
        self.opt_idx = idx
        print(f"costmap_cost: {costmap_costs[idx]}")
        print(f"goalcost: {goal_costs[idx]}")

        return self.candidate_arcs[idx], costs[idx], self.omegas[idx]
