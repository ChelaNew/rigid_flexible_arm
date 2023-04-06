from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

import simulation as s
import robot




def main():

	sim = s.simulation()
	anim = animation.FuncAnimation(sim.fig, sim.update, frames=1000, interval=20, blit=False)
	plt.show()


	# r = robot.flexible_link()
	# print(r.alphas)
	# print(r.evaluate_mode_functions(5))


if __name__ == "__main__":
	main()