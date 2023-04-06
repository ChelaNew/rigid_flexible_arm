from matplotlib import pyplot as plt
import numpy as np

import robot



class simulation:

	def __init__(self,num_links=5):

	fig = plt.figure()
	self.ax = fig.add_subplot(projection='3d')
	self.init_graph()

	self.lines = []
	for i in range(0,8):
		line, = self.ax.plot([], [], [], lw=2)
		self.lines.append(line)


	def init_graph(self):
		self.ax.axis("equal")

		self.ax.set_xlim3d([-10.0, 10.0])
		self.ax.set_xlabel('X')

		self.ax.set_ylim3d([-10.0, 10.0])
		self.ax.set_ylabel('Y')

		self.ax.set_zlim3d([0.0, 10.0])
		self.ax.set_zlabel('Z')

	def plot_len(line,rho,theta,phi,x_shift=0,y_shift=0,z_shift=0,plot_data=True):

		x = rho*np.cos(theta)*np.sin(phi) + x_shift
		y = rho*np.sin(theta)*np.sin(phi) + y_shift
		z = rho*np.cos(phi) + z_shift

		if (plot_data):
			line.set_data(x,y)
			line.set_3d_properties(z)

		return x[-1],y[-1],z[-1]

	


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# lines, = np.zeros()

# for line in lines:
# 	line = ax.plot(0,0,0)



def main():
	# r = robot.flexible_link()
	# print(r.alphas)
	# print(r.evaluate_mode_functions(5))

	r = rigid_links()

if __name__ == "__main__":
    main()