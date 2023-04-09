from matplotlib import pyplot as plt
import numpy as np

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


	def init_graph(self):
		self.ax.axis("equal")
		# self.ax.view_init(elev=90, azim=1)

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

	#position array has shape (num_joints, 4), or (num_joints, (x,y,z,1))
	def plot_rigid(self,position_arr):


		for i in range(len(position_arr)-1):
			x = np.linspace(position_arr[i,0],position_arr[i+1,0])
			y = np.linspace(position_arr[i,1],position_arr[i+1,1])
			z = np.linspace(position_arr[i,2],position_arr[i+1,2])

			self.lines[i].set_data(x,y)
			self.lines[i].set_3d_properties(z)



	def update(self,n):

		print(n)

		torque = 0.1
		# print(torque)

		rigi_pos_arr = self.rigid_links.run_dynamics(n,[0,0])
		self.plot_rigid(rigi_pos_arr)



		# rho_arr, theta_arr = self.FL_in_plane.run_dynamics_loop(n,0)
		# rho_arr, phi_arr = self.FL_out_plane.run_dynamics_loop(n,torque)

		# phi_arr = phi_arr + np.pi/2*np.ones(len(rho_arr))

		# self.plot_flexible(self.lines[0],rho_arr,theta_arr,phi_arr)

		# self.FL_in_plane.plots(['gen_coors','F'],n,i_stop=100)

			


















