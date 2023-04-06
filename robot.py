import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad
import numdifftools as nd


class flexible_link:

	def __init__(self):

		self.L = 5
		self.rho = 5.9292
		self.S = 0.06*0.04
		self.J = 5.12

		self.M_tot = self.rho*self.L*self.S
		self.M_bar = self.L/self.M_tot
		self.M_L = 0

		self.E = 5#6.21 * 10**7
		self.I = 1#??

		self.g = 9.18

		self.num_modes = 2
		self.alphas = self.find_eigenvalues()

		self.psi_L = self.evaluate_mode_functions(self.L)

		self.M = self.create_M(np.zeros(self.num_modes))
		self.C = self.create_C(0,np.zeros(self.num_modes),np.zeros(self.num_modes))
		self.K = self.create_K(np.zeros(self.num_modes))
		self.G = self.create_G(0)
		self.F = self.create_F(0,np.zeros(self.num_modes))
		self.Q = self.create_Q(0)


	def evaluate_mode_functions(self,x):

		psi_i = lambda x: (np.cosh(self.alphas*x) - np.cos(self.alphas*x)) - ((np.cosh(self.alphas*self.L) + np.cos(self.alphas*self.L)) / \
			(np.sinh(self.alphas*self.L) + np.sin(self.alphas*self.L))) * (np.sinh(self.alphas*x) - np.sin(self.alphas*x))

		return np.array(psi_i(x))



	def find_eigenvalues(self):

		#equation to deremine modal eigenvalues
		eigenvalues = lambda alpha_i: 1 + np.cos(alpha_i*self.L)*np.cosh(alpha_i*self.L) - \
			(self.M_L/(self.rho*self.S*self.L))*(alpha_i*self.L)*(np.sin(alpha_i*self.L)* \
			np.cosh(alpha_i*self.L) - np.cos(alpha_i*self.L)*np.sinh(alpha_i*self.L))

		#input range of possible eigenvalues into equation
		alpha = np.arange(0,5*self.num_modes,0.01)
		f = eigenvalues(alpha)

		#determine f values where the function changes sign for  numerical method
		sign_changes = np.where(np.diff(np.sign(f)))[0]

		#find corresponding x values
		sign_change_vals = []
		for i in sign_changes:
			sign_change_vals.append(alpha[i])

		#numerically find roots of equation to find eigenvalue values
		eigenvalues_array = []
		for i in range(len(sign_change_vals)-1):
			alpha_i = root_scalar(eigenvalues,bracket=[sign_change_vals[i], sign_change_vals[i+1]])
			eigenvalues_array.append(alpha_i.root)

		return np.array(eigenvalues_array[0:self.num_modes])


	def create_M(self,q):
		M = np.zeros((self.num_modes+1,self.num_modes+1))
		gammas = np.zeros(self.num_modes)

		for i in range(self.num_modes):
			gammas[i] = self.M_bar*quad(lambda x: x*self.evaluate_mode_functions(x)[i], 0, self.L)[0] + self.M_L*self.L*self.psi_L[i]
			M[i+1][i+1] = self.M_bar*quad(lambda x: self.evaluate_mode_functions(x)[i]**2, 0, self.L)[0] + self.M_tot*self.psi_L[i]**2

		m12 = gammas
		m21 = np.transpose(gammas)
		M[0,1:] = m12
		M[1:,0] = m21

		m11_sum1 = 0
		m11_sum2 = 0

		for i in range(self.num_modes):
			m11_sum1 += self.psi_L[i]**2 * q[i]**2
			m11_sum2 += quad(lambda x: self.evaluate_mode_functions(x)[i]**2 * q[i]**2, 0, self.L)[0]

		m11 = self.J + self.M_L*(self.L**2 + m11_sum1) + self.M_bar*m11_sum2
		M[0][0] = m11
			
		self.M = M


	def create_C(self,theta_dot,q,q_dot):
		C = np.zeros(self.num_modes+1)
	
		C0_sum1 = 0
		C0_sum2 = 0
		for i in range(self.num_modes):

			C0_sum1 += self.M_bar*q_dot[i]*q[i]*quad(lambda x: self.evaluate_mode_functions(x)[i]**2, 0, self.L)[0]
			C0_sum2 += self.M_tot*q_dot[i]*q[i]*self.psi_L[i]**2

			C[i+1] = -self.M_bar*quad(lambda x: self.evaluate_mode_functions(x)[i]**2 * theta_dot**2 * q[i] ,0, self.L)[0] + \
				self.M_tot*self.psi_L[i]**2 * theta_dot**2 * q[i]

		C[0] = 2*theta_dot*C0_sum1 + 2*theta_dot*C0_sum2

		self.C = C


	def create_K(self,q):
		K = np.zeros(self.num_modes+1)

		#in the paper this is a diagonal nxn matrix, but matrix eq solver needs it as a 1xn vector
		for i in range(self.num_modes):
			K[i+1] = self.E*self.I*quad(lambda x: nd.Derivative(self.evaluate_mode_functions, x, n=2)(x)[i], 0, self.L)[0]*q[i]

		self.K = K

	def create_G(self,theta):
		G = np.zeros(self.num_modes+1)

		G[0] = 0.5*self.M_tot*self.g*self.L*np.cos(theta) + self.M_L*self.g*self.L*np.cos(theta)

		self.G = G

	def create_F(self,theta_dot,q_dot,u=100,u1=0):
		F = np.zeros(self.num_modes+1)

		F[0] = u*theta_dot

		for i in range(self.num_modes):
			F[i+1] = 0.00*q_dot[i]

		F[1]=u1*q_dot[i]

		self.F = F

	def create_Q(self,torque_in):
		Q = np.transpose(np.zeros(self.num_modes+1))

		Q[0] = torque_in

		self.Q = Q


	

