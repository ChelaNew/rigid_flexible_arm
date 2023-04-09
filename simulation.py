from matplotlib import pyplot as plt
import numpy as np
import math

import flexible
import rigid



class simulation:

	def __init__(self,num_links=5):

		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(projection='3d')
		self.init_graph()

		self.lines = []
		self.line_index = 0
		for i in range(0,8):
			line, = self.ax.plot([], [], [], lw=2)
			self.lines.append(line)

		self.num_links = num_links

		self.rigid_links = rigid.rigid_links()
		self.FL_in_plane = flexible.flexible_link_dynamics()
		self.FL_out_plane = flexible.flexible_link_dynamics()
		self.FL_in_plane.set_initial_position(np.pi/2)

		self.num_joints = self.rigid_links.num_joints

		self.run_plots = True

		self.dt=1


	def init_graph(self):
		self.ax.axis("equal")
		self.ax.view_init(elev=90, azim=1)

		self.ax.set_xlim3d([-10.0, 10.0])
		self.ax.set_xlabel('X')

		self.ax.set_ylim3d([-10.0, 10.0])
		self.ax.set_ylabel('Y')

		self.ax.set_zlim3d([0.0, 10.0])
		self.ax.set_zlabel('Z')



	def plot_flexible(self,line,rho,theta,phi,x_shift=0,y_shift=0,z_shift=0,plot_data=True):

		x = rho*np.cos(theta)*np.sin(phi) + x_shift
		y = rho*np.sin(theta)*np.sin(phi) + y_shift
		z = rho*np.cos(phi) + z_shift


		if (plot_data):
			line.set_data(x,y)
			line.set_3d_properties(z)

		return x[-1],y[-1],z[-1]

	def convert_cart_to_sphere_coors(self,coors):
	

		r = np.sqrt(coors[0]**2 + coors[1]**2 + coors[2]**2)
		theta = np.arccos(coors[2]/r)
		phi = np.arctan2(coors[1],coors[0])

		return [r,theta,phi]

	def plot_line_cart(self,line,start_coor, end_coor):

		x = np.linspace(start_coor[0],end_coor[0])
		y = np.linspace(start_coor[1],end_coor[1])
		z = np.linspace(start_coor[2],end_coor[2])

		line.set_data(x,y)
		line.set_3d_properties(z)

	def plot_line_sphere(self,line,start_coor,end_coor):
		x_start = start_coor[0]*np.sin(start_coor[1])*np.cos(start_coor[2])
		y_start = start_coor[0]*np.sin(start_coor[1])*np.sin(start_coor[2])
		z_start = start_coor[0]*np.cos(start_coor[1])

		x_end = end_coor[0]*np.sin(end_coor[1])*np.cos(end_coor[2])
		y_end = end_coor[0]*np.sin(end_coor[1])*np.sin(end_coor[2])
		z_end = end_coor[0]*np.cos(end_coor[1])

		x = np.linspace(x_start,x_end)
		y = np.linspace(y_start,y_end)
		z = np.linspace(z_start,z_end)

		line.set_data(x,y)
		line.set_3d_properties(z)


	#position array has shape (num_joints, 4), or (num_joints, (x,y,z,1))
	def plot_rigid(self,position_arr):

		for i in range(len(position_arr)-1):
			x = np.linspace(position_arr[i,0],position_arr[i+1,0])
			y = np.linspace(position_arr[i,1],position_arr[i+1,1])
			z = np.linspace(position_arr[i,2],position_arr[i+1,2])

			self.lines[i].set_data(x,y)
			self.lines[i].set_3d_properties(z)

		return position_arr[-1,0],position_arr[-1,1],position_arr[-1,2],



	def update(self,n):

		torque = 0.1
		# print(torque)

		rigid_pos_arr, ee_accel = self.rigid_links.run_dynamics(n,[0,0.01],dt=self.dt)

		ee_accel = self.convert_cart_to_sphere_coors(ee_accel)

		rho_arr, theta_arr = self.FL_in_plane.run_dynamics_loop_single_link(n,0)
		rho_arr, phi_arr = self.FL_out_plane.run_dynamics_loop_single_link(n,0)


		x_end, y_end, z_end = self.plot_rigid(rigid_pos_arr)
		self.plot_flexible(self.lines[-1],rho_arr,theta_arr,phi_arr, x_shift=x_end, y_shift=y_end, z_shift=z_end)

		# self.FL_in_plane.plots(['gen_coors','F'],n,i_stop=100)

			


















