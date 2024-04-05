#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation  
import do_mpc
import casadi
from casadi import *
from scipy.interpolate import CubicSpline

class Point:
    """Helper class"""
    def __init__(self, x, y, s=0):
        self.x = x
        self.y = y
        self.s = s

class Path:
    """Helper class"""
    def __init__(self, x_pts, y_pts, s_pts, n_pts):
        assert len(x_pts) == len(y_pts)
        assert len(x_pts) == len(s_pts)
        self.x_pts = x_pts
        self.y_pts = y_pts
        self.s_pts = s_pts
        self.n_pts = n_pts
        self.points = [Point(x_pts[i], y_pts[i], s_pts[i]) for i in range(len(x_pts))]

class TrajectoryTransformer:
    
    def __init__(self, a=10.0, dt=0.01):
        self.a = a
        self.dt = dt
        self.path = self.generate_fig_eight(self.a)
        self.s = []
        self.d = []

        self.s_obs = []
        self.d_obs = []
        
        for point in self.path.points:
            s, d = self.frenet_transform(self.path, point.x, point.y)
            self.s.append(s)
            self.d.append(d)

        # traj_s_pts = np.arange(0.7, 30, 0.01)
        # traj_d_pts = 0.5 * np.power(np.e, -traj_s_pts/10)
        # self.x_traj = []
        # self.y_traj = []
        # # self.trajectory_path = self.generate_exponential_trajectory(1.0, 1.0)
        # for index in range(len(traj_s_pts)):
        #     x, y = self.inverse_frenet_transform(traj_s_pts[index], traj_d_pts[index])
        #     self.x_traj.append(x)
        #     self.y_traj.append(y)

        """
        obs_x, obs_y = self.generate_circular_obstacles(x_c=self.a*0.6, y_c=0.0, radius=self.a*0.3)
        self.obstacles.add_obstacles(obs_x, obs_y)

        obs_x, obs_y = self.generate_circular_obstacles(x_c=-self.a*0.6, y_c=0.0, radius=self.a*0.3)
        self.obstacles.add_obstacles(obs_x, obs_y)

        for obstacle_point in self.obstacles.obstacle_points:
            s_obs, d_obs = self.frenet_transform(self.path, obstacle_point.x, obstacle_point.y)
            self.s_obs.append(s_obs)
            self.d_obs.append(d_obs)

        fig, (ax_c, ax_f) = plt.subplots(1, 2)
        marker_size = 5
        self.c = []
        # self.c = [[0 / len(self.path.x_pts), 0 / len(self.path.x_pts), i / len(self.path.x_pts)] for i in range(len(self.path.x_pts))]
        for i in range(len(self.path.x_pts)):
            if np.tan(self.path.n_pts[i]) > 0:
                self.c.append([0, 1, 0])
            else:
                self.c.append([1, 0, 0])
        ax_c.scatter(self.path.x_pts, self.path.y_pts, color=self.c, s=marker_size, label='reference')
        ax_c.scatter(self.x_traj, self.y_traj, color='orange', s=marker_size, label='traj inverse')
        # ax_c.scatter(self.obstacles.x_pts, self.obstacles.y_pts, color='r', s=marker_size, label='obstacle')
        ax_f.scatter(traj_s_pts, traj_d_pts, color='g', s=marker_size, label='traj')
        ax_f.scatter(self.s, self.d, color=self.c, s=marker_size, label='reference')
        ax_f.scatter(self.s_obs, self.d_obs, color='r', s=marker_size, label='obstacle')

        reverse_x = []
        reverse_y = []
        for i in range(len(self.s)):
            x, y = self.inverse_frenet_transform(self.s[i], self.d[i])
            reverse_x.append(x)
            reverse_y.append(y)
        # ax_c.plot(reverse_x, reverse_y, label='inverse', alpha=0.5, color='y')
        # ax_c.plot(reverse_x, reverse_y, s=marker_size * 1.5, label='inverse', alpha=0.5, color='y')
        
        ax_c.legend()
        ax_f.legend()
        plt.show()
        """
        
    def generate_fig_eight(self, a):
        """Generates points, progress variable s and slope for figure 8 trajectory"""
        t = np.arange(0, 2 * np.pi + self.dt, self.dt)
        x = a * np.sin(t)
        y = a * np.sin(t) * np.cos(t)
        s = a * np.sqrt(np.cos(t)**2 + np.cos(2*t)**2) * self.dt
        n = np.arctan2(np.cos(2*t), np.cos(t))
        for i in range(1, len(s)):
            s[i] = s[i] + s[i-1]
        return Path(x, y, s, n)

    def get_distance(self, point_a: Point, point_b: Point):
        """Returns Eucledian distance"""
        return np.sqrt(
            (point_a.x - point_b.x) ** 2 + (point_a.y - point_b.y)**2
        )
    
    def get_closest_point(self, path: Path, ref_point):
        """Helper function for the frenet transform"""
        closest_distance = 10000
        closest_index = -1
        for index, point in enumerate(path.points):
            distance = self.get_distance(point, ref_point)
            if distance < closest_distance:
                closest_distance = distance
                closest_index = index

        return (path.points[closest_index], closest_index)
    
    def get_inclination(self, point_a: Point, point_b: Point):
        """Returns the slope of the line joining two input points"""
        return np.arctan2((point_a.y - point_b.y), (point_a.x - point_b.x))
    
    def frenet_transform(self, path: Path, input_x, input_y):
        """Transforms cartesian coordinates to frenet coordinates"""
        s = 0
        d = 0
        if path == []:
            path = self.path
        closest_point, closest_point_index = self.get_closest_point(path, Point(input_x, input_y, 0.0))
        s = path.s_pts[closest_point_index]
        d = self.get_distance(closest_point, Point(input_x, input_y, 0.0))
        theta = path.n_pts[closest_point_index]
        psi = self.get_inclination(closest_point, Point(input_x, input_y, 0.0))
        d_sign = -1 if psi < theta else +1
        d = d_sign * d
        return (s, d)
    
    def inverse_frenet_transform(self, input_s, input_d):
        """Transforms frenet coordinates to cartesian coordinates"""
        min_dist = 1000
        min_index = -1
        for index, s_pt in enumerate(self.path.s_pts):
            dist = abs(input_s - s_pt)
            if dist < min_dist:
                min_index = index
                min_dist = dist
        intermediate_point_index = min_index
        intermediate_point = self.path.points[intermediate_point_index]
        inclination = self.path.n_pts[intermediate_point_index]
        inclination_sign = +1 if input_d > 0 else -1
        # ninety_plus_theta = np.arctan2(path.n_pts[intermediate_point_index], 1)
        ninety_plus_theta = inclination + inclination_sign * np.pi/2
        x = intermediate_point.x + abs(input_d) * np.cos(ninety_plus_theta)
        y = intermediate_point.y + abs(input_d) * np.sin(ninety_plus_theta)
        return (x, y)

class RobotSimulator:

    def __init__(self, dt=0.1, a=5.0):

        self.a = a # Scaling parameter of figure 8 trajectory
        self.trajectory_transformer = TrajectoryTransformer(a=self.a, dt=0.01)
        self.ddot_max = 0.5
        self.desired_velocity = 1.0

        # Defining model
        self.dt = dt
        self.mpc_model_type = 'discrete'
        self.mpc_model = do_mpc.model.Model(self.mpc_model_type)

        # Internal model states
        state_s = self.mpc_model.set_variable(var_name='state_s', var_type='_x')
        state_d = self.mpc_model.set_variable(var_name='state_d', var_type='_x')
        state_s_dot = self.mpc_model.set_variable(var_name='state_s_dot', var_type='_x')
        state_d_dot = self.mpc_model.set_variable(var_name='state_d_dot', var_type='_x')

        input_s_ddot = self.mpc_model.set_variable(var_name='input_s_ddot', var_type='_u')        
        input_d_ddot = self.mpc_model.set_variable(var_name='input_d_ddot', var_type='_u')

        state_s_next = state_s + state_s_dot * dt + 0.5 * input_s_ddot * dt**2
        state_d_next = state_d + state_d_dot * dt + 0.5 * input_d_ddot * dt**2
        state_s_dot_next = state_s_dot + input_s_ddot * dt
        state_d_dot_next = state_d_dot + input_d_ddot * dt

        self.mpc_model.set_rhs('state_s', state_s_next)
        self.mpc_model.set_rhs('state_d', state_d_next)
        self.mpc_model.set_rhs('state_s_dot', state_s_dot_next)
        self.mpc_model.set_rhs('state_d_dot', state_d_dot_next)

        # Cost function
        deviation_cost = 15 * state_d ** 2
        progress_cost = - 5 * state_s
        desired_velocity_cost = (state_s_dot - self.desired_velocity)**2
        # desired_velocity_cost = SX(0.0)
        self.mpc_model.set_expression('deviation_cost', deviation_cost)
        self.mpc_model.set_expression('progress_cost', progress_cost)
        self.mpc_model.set_expression('desired_velocity_cost', desired_velocity_cost)
        self.mpc_model.setup()

        self.n_horizon = 5 * int(1 / self.dt)

        self.mpc = do_mpc.controller.MPC(self.mpc_model)
        solver = 'ma97'
        setup_mpc = {
            'n_horizon': self.n_horizon,
            't_step': self.dt,
            'store_full_solution': True,
            'nlpsol_opts': {'ipopt.linear_solver': solver,'ipopt.print_level': 0}
        }
        self.mpc.set_param(**setup_mpc)
        lterm = self.mpc_model.aux['deviation_cost'] + self.mpc_model.aux['progress_cost'] + self.mpc_model.aux['desired_velocity_cost']
        mterm = self.mpc_model.aux['deviation_cost'] + self.mpc_model.aux['progress_cost'] + self.mpc_model.aux['desired_velocity_cost']
        self.mpc.set_objective(lterm=lterm, mterm=mterm)
        rterm = 1e-4
        self.mpc.set_rterm(input_s_ddot=rterm, input_d_ddot=rterm)

        # State Constraints
        self.mpc.bounds['lower', '_x', 'state_s_dot'] = -self.desired_velocity
        self.mpc.bounds['upper', '_x', 'state_s_dot'] = self.desired_velocity
        # Input Constraints
        self.mpc.bounds['lower', '_u', 'input_s_ddot'] = -self.ddot_max
        self.mpc.bounds['upper', '_u', 'input_s_ddot'] = self.ddot_max
        self.mpc.bounds['lower', '_u', 'input_d_ddot'] = -self.ddot_max
        self.mpc.bounds['upper', '_u', 'input_d_ddot'] = self.ddot_max

        self.obstacle_s = 10.0
        self.obstacle_d = 0.0
        self.obstacle_radius = 0.5
        # Obstacle constraint
        # self.mpc.set_nl_cons('obstacle_cons', -((state_s - self.obstacle_s) ** 2 + (state_d - self.obstacle_d)**2), -(self.obstacle_radius ** 2))

        self.mpc.setup()
        self.mpc.x0 = np.array([0.0, 0.0, 0.0, 0.0])  # s, d, s_dot, d_dot
        self.mpc.u0 = np.array([0.0, 0.0]) # s_ddot, d_ddot

        self.mpc.set_initial_guess()

    def visualize(self, n_steps):
        marker_size = 1
        
        fig = plt.figure(figsize=(16, 9))
        ax_c = fig.add_subplot(3, 5, (1, 13))
        ax_f = fig.add_subplot(3, 5, (4, 5))
        ax_s = fig.add_subplot(3, 5, (9, 10))
        ax_d = fig.add_subplot(3, 5, (14, 15))

        time = np.arange(0, (self.n_horizon + 1) * self.dt, self.dt)

        self.anim_traj_s = ax_s.plot(time, self.s_plot[0], label="$s(t)$")[0]
        self.anim_traj_s_dot = ax_s.plot(time, self.s_dot_plot[0], label="$\dot s(t)$")[0]
        self.anim_traj_s_ddot = ax_s.plot(time, self.s_ddot_plot[0], label="$\ddot s(t)$")[0]

        self.anim_traj_d = ax_d.plot(time, self.d_plot[0], label="$d(t)$")[0]
        self.anim_traj_d_dot = ax_d.plot(time, self.d_dot_plot[0], label="$\dot d(t)$")[0]
        self.anim_traj_d_ddot = ax_d.plot(time, self.d_ddot_plot[0], label="$\ddot d(t)$")[0]

        # obs_x, obs_y = self.trajectory_transformer.inverse_frenet_transform(self.obstacle_s, self.obstacle_d)
        # obstacle_circle = plt.Circle((obs_x, obs_y), self.obstacle_radius, color='r', fill=False)
        # ax_c.add_patch(obstacle_circle)       # static
        
        self.anim_c_robot = ax_c.scatter([0.0], [0.0], color='blue', s=75, label='robot')
        self.anim_c_traj = ax_c.plot(self.x_plot[0], self.y_plot[0], color='lawngreen', label='optimal trajectory', linewidth=3)[0]
        self.anim_c_ref = ax_c.scatter(self.trajectory_transformer.path.x_pts, self.trajectory_transformer.path.y_pts, s=marker_size, color='green', alpha=0.75, label='reference') # static

        # self.anim_f_obs = ax_f.scatter(self.obstacle_s, self.obstacle_d, color='red', s=100, label='obstacle')  # static
        self.anim_f_robot = ax_f.scatter([0.0], [0.0], color='blue', label='robot', s=75)
        self.anim_f_traj = ax_f.plot(self.s_plot[0], self.d_plot[0], color='lawngreen', label='optimal trajectory', linewidth=3)[0]
        self.anim_f_ref = ax_f.scatter(self.trajectory_transformer.s, self.trajectory_transformer.d, s=marker_size, color='green', alpha=0.75, label='reference') # static

        ax_c.set_xlabel('Position x [m]')
        ax_c.set_ylabel('Position y [m]')

        ax_f.set_xlabel('Progress s [m]')
        ax_f.set_ylabel('Deviation d [m]')

        ax_s.set_xlabel('Time [s]')
        ax_s.set_ylabel('Progress s [m]')

        ax_d.set_xlabel('Time t [s]')
        ax_d.set_ylabel('Deviation d [m]')

        ax_c.set_xlim(left=-1.25 * self.a , right=1.25 * self.a)
        ax_c.set_ylim(bottom=-1.25 * self.a, top=1.25 * self.a)

        ax_f.set_ylim(bottom=-5, top=5)
        ax_f.set_xlim(left=-1, right=35)

        ax_c.set_aspect('equal')
        ax_f.set_aspect('equal')

        ax_c.legend()
        ax_f.legend(loc='upper right')
        ax_s.legend()
        ax_d.legend()
        # ax_s.set_xlim(left=0, right=60)
        ax_s.set_ylim(bottom=0, top=35)

        anim = animation.FuncAnimation(fig, self.anim_update, frames=n_steps, interval=100, blit=True)
        writer = animation.PillowWriter(fps=10, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('figure_eight.gif', writer)
        plt.show()
    
    def simulate(self, n_steps):

        fig = plt.figure(figsize=(16, 9))
        ax_c = fig.add_subplot(3, 5, (1, 13))
        ax_f = fig.add_subplot(3, 5, (4, 5))
        # ax_curv = fig.add_subplot(3, 5, (9, 10))
        ax_s = fig.add_subplot(3, 5, (9, 10))
        ax_d = fig.add_subplot(3, 5, (14, 15))
        marker_size = 1
        
        # Initialize lists to store simulation data
        self.s_plot = []
        self.s_dot_plot = []
        self.s_ddot_plot = []
        self.d_plot = []
        self.d_dot_plot = []
        self.d_ddot_plot = []

        self.x_plot = []
        self.y_plot = []
        self.theta_plot = []

        
        x0 = np.array([2.0, 3.0, 0.0, 0.0]) # Initial condition in s and d
        for i in range(n_steps):

            u0 = self.mpc.make_step(x0)
            # mpc gives trajectory in s, d
            s_traj = self.mpc.data.prediction(('_x', 'state_s')).flatten()
            d_traj = self.mpc.data.prediction(('_x', 'state_d')).flatten()
            s_dot_traj = self.mpc.data.prediction(('_x', 'state_s_dot')).flatten()
            d_dot_traj = self.mpc.data.prediction(('_x', 'state_d_dot')).flatten()
            s_ddot_traj = np.append(np.array([0.0]), self.mpc.data.prediction(('_u', 'input_s_ddot')).flatten())
            d_ddot_traj = np.append(np.array([0.0]), self.mpc.data.prediction(('_u', 'input_d_ddot')).flatten())
    

            x_pts = []
            y_pts = []
            theta_pts = [0.0]
            x0 = np.array([s_traj[1], d_traj[1], s_dot_traj[1], d_dot_traj[1]])
            # print(f"shape: {s_traj.shape[0]}")

            for i in range(s_traj.shape[0]):
                x, y = self.trajectory_transformer.inverse_frenet_transform(s_traj[i], d_traj[i])
                x_pts.append(x)
                y_pts.append(y)
                if(i > 0):
                    theta_pts.append(np.arctan2(y_pts[-1] - y_pts[-2], x_pts[-1] - x_pts[-2]))

            # print(self.n_horizon, self.dt)
            # curvature = self.check_path_feasibility(x_pts, y_pts, self.n_horizon, self.dt)
            # print("Min turning radius: ", min(curvature[:5]))

            # Store all values
            self.s_plot.append(s_traj)
            self.s_dot_plot.append(s_dot_traj)
            self.s_ddot_plot.append(s_ddot_traj)

            self.d_plot.append(d_traj)
            self.d_dot_plot.append(d_dot_traj)
            self.d_ddot_plot.append(d_ddot_traj)

            self.x_plot.append(x_pts)
            self.y_plot.append(y_pts)
            self.theta_plot.append(theta_pts)

            """
            # Debugging
            ax_c.clear()
            ax_f.clear()
            ax_curv.clear()
            ax_d.clear()
            ax_s.clear()
            
            time = np.arange(0, (self.n_horizon + 1) * self.dt, self.dt)
            ax_curv.plot(time, curvature)

            ax_s.plot(time, s_traj, label="$s(t)$")
            ax_s.plot(time, s_dot_traj, label="$\dot s(t)$")
            ax_s.plot(time, s_ddot_traj, label="$\ddot s(t)$")

            ax_d.plot(time, d_traj, label="$d(t)$")
            ax_d.plot(time, d_dot_traj, label="$\dot d(t)$")
            ax_d.plot(time, d_ddot_traj, label="$\ddot d(t)$")

            obs_x, obs_y = self.trajectory_transformer.inverse_frenet_transform(self.obstacle_s, self.obstacle_d)
            obstacle_circle = plt.Circle((obs_x, obs_y), self.obstacle_radius, color='r', fill=False)
            # ax_c.scatter(obs_x, obs_y, color='red', s=100, label='obstacle')
            ax_c.add_patch(obstacle_circle)        
            ax_c.scatter(x_pts[0], y_pts[0], color='black', s=75, label='robot')

            ax_c.plot(x_pts, y_pts, color='lawngreen', label='traj', linewidth=3)
            ax_c.scatter(self.trajectory_transformer.path.x_pts, self.trajectory_transformer.path.y_pts, s=marker_size, color='black', alpha=0.75, label='ref')

            ax_f.scatter(self.obstacle_s, self.obstacle_d, color='red', s=100, label='obstacle')
            ax_f.scatter(s_traj[0], d_traj[0], color='black', label='robot', s=75)
            ax_f.plot(s_traj, d_traj, color='lawngreen', label='traj', linewidth=3)
            ax_f.scatter(self.trajectory_transformer.s, self.trajectory_transformer.d, s=marker_size, color='black', alpha=0.75, label='ref')
            
            ax_c.legend()
            ax_f.legend()
            ax_s.legend()
            ax_d.legend()

            ax_c.set_aspect('equal')
            ax_f.set_aspect('equal')

            ax_c.set_xlim(left=-self.a, right=self.a)
            ax_c.set_ylim(bottom=-self.a, top=self.a)
            ax_f.set_ylim(bottom=-4, top=4)
            
            plt.draw()
            plt.pause(0.001)
            """
                       

        print(f"Simulation for {n_steps} complete")

    def anim_update(self, timestep):
        time = np.arange(0, (self.n_horizon + 1) * self.dt, self.dt)

        self.anim_traj_s.set_ydata(self.s_plot[timestep])
        self.anim_traj_s_dot.set_ydata(self.s_dot_plot[timestep])
        self.anim_traj_s_ddot.set_ydata(self.s_ddot_plot[timestep])
        self.anim_traj_s.set_xdata(time)
        self.anim_traj_s_dot.set_xdata(time)
        self.anim_traj_s_ddot.set_xdata(time)

        self.anim_traj_d.set_ydata(self.d_plot[timestep])
        self.anim_traj_d_dot.set_ydata(self.d_dot_plot[timestep])
        self.anim_traj_d_ddot.set_ydata(self.d_ddot_plot[timestep])
        self.anim_traj_d.set_xdata(time)
        self.anim_traj_d_dot.set_xdata(time)
        self.anim_traj_d_ddot.set_xdata(time)

        self.anim_c_robot.set_offsets([self.x_plot[timestep][0], self.y_plot[timestep][0]])
        self.anim_c_traj.set_xdata(self.x_plot[timestep])
        self.anim_c_traj.set_ydata(self.y_plot[timestep])

        self.anim_f_robot.set_offsets([self.s_plot[timestep][0], self.d_plot[timestep][0]])
        self.anim_f_traj.set_xdata(self.s_plot[timestep])
        self.anim_f_traj.set_ydata(self.d_plot[timestep])

        return self.anim_traj_s, self.anim_traj_s_dot, self.anim_traj_s_ddot, self.anim_traj_d, self.anim_traj_d_dot, self.anim_traj_d_ddot, self.anim_c_robot, self.anim_c_traj, self.anim_c_ref, self.anim_f_robot, self.anim_f_ref, self.anim_f_traj, 
        # return self.anim_traj_s, self.anim_traj_s_dot, self.anim_traj_s_ddot, self.anim_traj_d, self.anim_traj_d_dot, self.anim_traj_d_ddot, self.anim_c_robot, self.anim_c_traj, self.anim_c_ref, self.anim_f_obs, self.anim_f_robot, self.anim_f_ref, self.anim_f_traj, 

    def check_path_feasibility(self, x_pts, y_pts, n_horizon, dt):
        """Checks for curvature of the path"""
        curvature = []
        t_horizon = (n_horizon + 1) * dt  
        t = np.arange(0, t_horizon, dt)
        x_spline = CubicSpline(t, x_pts)
        y_spline = CubicSpline(t, y_pts)
        x_dot = x_spline.derivative(1)
        y_dot = y_spline.derivative(1)
        x_ddot = x_spline.derivative(2)
        y_ddot = y_spline.derivative(2)
        curvature = []
        for t_ in t:
            c = abs(x_ddot(t_) * y_dot(t_) - x_dot(t_) * y_ddot(t_)) / np.power((x_dot(t_) ** 2 + y_dot(t_) ** 2), 3/2)
            curvature.append(1/c)
        return curvature



if __name__ == '__main__':
    a = 1.0
    simulator = RobotSimulator()
    n_steps = 300
    simulator.simulate(n_steps=n_steps)
    simulator.visualize(n_steps=n_steps)
