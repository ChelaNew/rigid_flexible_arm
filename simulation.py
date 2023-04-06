from matplotlib import pyplot as plt
import numpy as np

import robot
import rigid


class simulation:

	def __init__(self,num_links=5):

		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(projection='3d')
		self.init_graph()

		self.lines = []
		for i in range(0,8):
			line, = self.ax.plot([], [], [], lw=2)
			self.lines.append(line)

		self.num_links = num_links

		self.rigid_links = rigid.rigid_links()
		self.num_joints = self.rigid_links.num_joints


	def init_graph(self):
		self.ax.axis("equal")

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

		#add in origin as starting point
		position_arr = np.concatenate([np.zeros((1,4)),position_arr],axis=0)


		for i in range(self.num_links):
			x = np.linspace(position_arr[i][0],position_arr[i+1][0])
			y = np.linspace(position_arr[i][1],position_arr[i+1][1])
			z = np.linspace(position_arr[i][2],position_arr[i+1][2])

			self.lines[i].set_data(x,y)
			self.lines[i].set_3d_properties(z)



	def update(self,n):
		torque = [100,0,0,0,0]

		self.rigid_links.run_dynamics(n,torque)
		self.rigid_links.convert_coordinate_frames(n)
		rigid_link_pos = self.rigid_links.convert_coordinate_frames(n)

		self.plot_rigid(rigid_link_pos)


















