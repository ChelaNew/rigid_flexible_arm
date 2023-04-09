from matplotlib import pyplot as plt
import numpy as np

from scipy.optimize import root_scalar
from scipy.integrate import quad
import numdifftools as nd
from scipy import linalg


class flexible_link_dynamics:

	def __init__(self,num_iter=1000):

		self.L = 5
		self.rho = 5.9292
		self.S = 0.06*0.04
		self.J = 5.12

		self.M_tot = self.rho*self.L*self.S
		self.M_bar = self.M_tot/self.L
		self.M_L = 0

		self.E = 5#6.21 * 10**7
		self.I = 0.0001#??

		self.g = 9.18

		self.num_modes = 1
		self.alphas = self.find_modes(self.num_modes)

		self.psi_L = self.psi(self.L)

		self.theta = np.zeros(num_iter)
		self.theta_dot = np.zeros(num_iter)
		self.theta_ddot = np.zeros(num_iter)

		self.q = np.zeros((self.num_modes,num_iter))
		self.q_dot = np.zeros((self.num_modes,num_iter))
		self.q_ddot = np.zeros((self.num_modes,num_iter))

		self.torque_in = 0

		self.M = np.zeros((self.num_modes+1,self.num_modes+1))
		self.Q = np.zeros(self.num_modes+1)
		self.C = np.zeros(self.num_modes+1)
		self.K = np.zeros(self.num_modes+1)
		self.G = np.zeros(self.num_modes+1)
		self.F = np.zeros(self.num_modes+1)

		self.M_arr = np.zeros((self.num_modes+1,self.num_modes+1,num_iter))
		self.Q_arr = np.zeros((self.num_modes+1,num_iter))
		self.C_arr = np.zeros((self.num_modes+1,num_iter))
		self.K_arr = np.zeros((self.num_modes+1,num_iter))
		self.G_arr = np.zeros((self.num_modes+1,num_iter))
		self.F_arr = np.zeros((self.num_modes+1,num_iter))


	def set_initial_position(self,theta):
		self.theta[0] = theta


	def psi(self,x):

		psi_i = (np.cosh(self.alphas*x) - np.cos(self.alphas*x)) - ((np.cosh(self.alphas*self.L) + np.cos(self.alphas*self.L)) \
			/ (np.sinh(self.alphas*self.L) + np.sin(self.alphas*self.L))) * \
			(np.sinh(self.alphas*x) - np.sin(self.alphas*x))

		return np.array(psi_i)


	def find_modes(self,n):

		f = lambda alpha_i: 1 + np.cos(alpha_i*self.L)*np.cosh(alpha_i*self.L) - \
		(self.M_L/(self.rho*self.S*self.L))*(alpha_i*self.L)*(np.sin(alpha_i*self.L)*np.cosh(alpha_i*self.L) - \
			np.cos(alpha_i*self.L)*np.sinh(alpha_i*self.L))

		x = np.arange(0,10,0.01)
		y = f(x)

		zero_crossings = np.where(np.diff(np.sign(y)))[0]

		vals = []
		for i in zero_crossings:
			vals.append(x[i])

		alphas = []
		for i in range(len(vals)-1):
			alpha_i = root_scalar(f,bracket=[vals[i], vals[i+1]])
			alphas.append(alpha_i.root)

		return np.array(alphas[0:n])


	def create_M(self,q):

		M = np.zeros((self.num_modes+1,self.num_modes+1))

		gammas = np.zeros(self.num_modes)
		for i in range(self.num_modes):
			gammas[i] = self.M_bar*quad(lambda x: x*self.psi(x)[i], 0, self.L)[0] + self.M_L*self.L*self.psi_L[i]
			M[i+1][i+1] = self.M_bar*quad(lambda x: self.psi(x)[i]**2, 0, self.L)[0] + self.M_bar*self.psi_L[i]**2

		m12 = gammas
		m21 = np.transpose(gammas)
		M[0,1:] = m12
		M[1:,0] = m21

		m11_sum1 = 0
		m11_sum2 = 0

		for i in range(self.num_modes):
			m11_sum1 += self.psi_L[i]**2 * q[i]**2
			m11_sum2 += quad(lambda x: self.psi(x)[i]**2 * q[i]**2, 0, self.L)[0]

		m11 = self.J + self.M_L*(self.L**2 + m11_sum1) + self.M_bar*m11_sum2

		M[0][0] = m11
		
		return M


	def create_C(theta_dot,q,q_dot,alphas):
		C = np.zeros(1+num_modes)
		
		C0_sum1 = 0
		C0_sum2 = 0
		for i in range(num_modes):

			C0_sum1 += m_bar*q_dot[i]*q[i]*quad(lambda x: psi(x,alphas)[i]**2, 0, L)[0]
			C0_sum2 += m*q_dot[i]*q[i]*psi_L[i]**2

			C[i+1] = -m_bar*quad(lambda x: psi(x,alphas)[i]**2 * theta_dot**2 * q[i] ,0, L)[0] + m*psi_L[i]**2 * theta_dot**2 * q[i]

		C[0] = 2*theta_dot*C0_sum1 + 2*theta_dot*C0_sum2

		return C


	def create_K(self,q):
		K = np.zeros(self.num_modes+1)

		for i in range(self.num_modes):

			K[i+1] = self.E*self.I*quad(lambda x: nd.Derivative(self.psi, x, n=2)(x)[i], 0, self.L)[0]*q[i]

		return K	

	def create_G(theta):
		G = np.zeros(num_modes+1)

		G[0] = 0.5*m*g*L*np.cos(theta) + m_p*g*L*np.cos(theta)

		return G

	def create_F(theta_dot,q_dot,u=100,u1=0):
		F = np.zeros(num_modes+1)

		F[0] = u*theta_dot

		for i in range(num_modes):
			F[i+1] = 0.00*q_dot[i]

		F[1]=u1*q_dot[i]

		return F


	def create_Q(self,theta_tau):

		Q = np.transpose(np.zeros(self.num_modes+1))

		Q[0] = theta_tau

		return Q

	def generate_plotting_arrays(self,i,dx=0.1):
		position_arr = np.arange(0,self.L,dx)
		angle_arr = [0]

		for x in np.arange(dx,self.L,dx):
			phi = self.psi(x)

			displacement = np.dot(phi,self.q[:,i])

			angle_arr.append( np.arctan(displacement/x) )

		angle_arr += self.theta[i]

		return position_arr,angle_arr


	def run_dynamics_loop_single_link(self,i,torque_in,dt=1):

		self.torque_in = torque_in

		self.M = self.create_M(self.q[:,i])
		self.Q = self.create_Q(torque_in)
		# self.C = self.create_C(self.theta_dot[i],self.q[:,i],self.q_dot[:,i])
		self.K = self.create_K(self.q[:,i])
		# self.G = self.create_G(self.theta[i])
		# self.F = self.create_F(self.theta_dot[i],self.q_dot[:,i])

		accel_arr = linalg.solve(self.M, self.Q-self.K)
		self.theta_ddot[i+1] = accel_arr[0]
		self.q_ddot[:,i+1] = accel_arr[1:]


		self.theta_dot[i+1] = self.theta_dot[i] + self.theta_ddot[i+1]*dt
		self.theta[i+1] = self.theta[i] + self.theta_dot[i+1]*dt

		self.q_dot[:,i+1] = self.q_dot[:,i] + self.q_ddot[:,i+1]*dt
		self.q[:,i+1] = self.q[:,i] + self.q_dot[:,i+1]*dt

		# print(self.theta[i+1],self.q[0][i+1])

		return self.generate_plotting_arrays(i)

	def run_dynamics_loop_with_rigid(self,i,torque_in,dt=1,external_accel=0):

		self.torque_in = torque_in

		self.M = self.create_M(self.q[:,i])
		self.Q = self.create_Q(torque_in)
		# self.C = self.create_C(self.theta_dot[i],self.q[:,i],self.q_dot[:,i])
		self.K = self.create_K(self.q[:,i])
		# self.G = self.create_G(self.theta[i])
		# self.F = self.create_F(self.theta_dot[i],self.q_dot[:,i])

		accel_arr = linalg.solve(self.M, self.Q-self.K)
		self.theta_ddot[i+1] = accel_arr[0] + external_accel
		self.q_ddot[:,i+1] = accel_arr[1:]


		self.theta_dot[i+1] = self.theta_dot[i] + self.theta_ddot[i+1]*dt
		self.theta[i+1] = self.theta[i] + self.theta_dot[i+1]*dt

		self.q_dot[:,i+1] = self.q_dot[:,i] + self.q_ddot[:,i+1]*dt
		self.q[:,i+1] = self.q[:,i] + self.q_dot[:,i+1]*dt

		# print(self.theta[i+1],self.q[0][i+1])

		return self.generate_plotting_arrays(i)




	def plots(self,type,i=0,i_stop=100,dx=0.1):

		if 'modes' in type:
			phi = []

			for x in np.arange(dx,self.L,dx):
				phi.append(self.evaluate_mode_functions(x))

			plt.plot(phi)

		if 'gen_coors' in type:

			if (i >= i_stop):
				plt.figure(2)

				plt.subplot(self.num_modes+1, 3, 1)
				plt.plot(self.theta[:i])
				plt.subplot(self.num_modes+1, 3, 2)
				plt.plot(self.theta_dot[:i])
				plt.subplot(self.num_modes+1, 3, 3)
				plt.plot(self.theta_ddot[:i])

				for j in range(1,self.num_modes+1):
					plt.subplot(self.num_modes+1, 3, 3*j+1)
					plt.plot(self.q[j-1][:i])
					plt.subplot(self.num_modes+1, 3, 3*j+2)
					plt.plot(self.q_dot[j-1][:i])
					plt.subplot(self.num_modes+1, 3, 3*j+3)
					plt.plot(self.q_ddot[j-1][:i])

		if 'M' in type:
			self.M_arr[:,:,i] = self.M

			if (i >= i_stop):
				plt.figure(3)

				plot_index=1
				for j in range(self.num_modes+1):
					for k in range(self.num_modes+1):
						plt.subplot(self.num_modes+1, self.num_modes+1, plot_index)
						plt.plot(self.M_arr[j][k][:i])
						plot_index += 1


		if 'F' in type:
			self.F_arr[:,i] = self.F

			if (i >= i_stop):
				plt.figure(4)

				for j in range(self.num_modes+1):
					plt.subplot(1,self.num_modes+1, j+1)
					plt.plot(self.F_arr[j][:i])

		if (i >= i_stop):
			plt.show()



			





	

