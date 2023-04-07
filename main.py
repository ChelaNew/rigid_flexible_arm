from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

import flexible as f
import simulation as s




def main():

	# flex = f.flexible_link()
	# flex.graph_modes()


	sim = s.simulation()
	anim = animation.FuncAnimation(sim.fig, sim.update, frames=1000, interval=20, blit=False)
	plt.show()


if __name__ == "__main__":
	main()