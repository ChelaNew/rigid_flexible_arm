import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import animation


class rigid_links:
	def __init__(self,num_iter=1000):
		self.num_joints = 5
		self.link_length = np.array([1.5,10.1,5.5,6.0,10.5])  # units cm
		self.link_radius = np.array([1,1,1,1,1])  # units cm
		self.link_mass = np.array([20,20,20,20,20])   # units grams
		self.s_hat = np.array([[0,0,.75],[5.05,0,0],[4.75,0,0],[0,3.0,0],[0,0,5.25]])

		self.DH_a = np.array([0,0,101,55,0])
		self.DH_alpha = np.array([0,-np.pi/2,0,0,np.pi/2])
		self.DH_d = np.array([65,0,0,0,60])
		self.DH_theta = np.array([0,-np.pi/2,0,np.pi/2,0])


		self.H_hat = np.zeros((self.num_joints,4,4))
		self.calc_H_hat()

		self.q = np.zeros((num_iter,self.num_joints),dtype=float)
		self.q_dot = np.zeros((num_iter,self.num_joints),dtype=float)
		self.q_ddot = np.zeros((num_iter,self.num_joints),dtype=float)
		


	def calc_Tk(self,k,theta):
	
		assert (k >= 1) and (k <= self.num_joints), "Error: calc_Tk - The k value is out of range."
		
		# Define the DH parameters (all indexed from zero)
		# Note that z-axes are reversed from above (theta and alpha are multiplied by -1)\

		DH_a = self.DH_a
		DH_alpha = self.DH_alpha
		DH_d = self.DH_d
		DH_theta = self.DH_theta

		
		Tk = np.matrix([
			[np.cos(DH_theta[k-1]+theta), -np.sin(DH_theta[k-1]+theta), 0, DH_a[k-1]],
			[np.cos(DH_alpha[k-1]) * np.sin(DH_theta[k-1]+theta), np.cos(DH_alpha[k-1]) * np.cos(DH_theta[k-1]+theta),-np.sin(DH_alpha[k-1]),-np.sin(DH_alpha[k-1]) * DH_d[k-1]],
			[np.sin(DH_alpha[k-1]) * np.sin(DH_theta[k-1]+theta), np.sin(DH_alpha[k-1]) * np.cos(DH_theta[k-1]+theta),np.cos(DH_alpha[k-1]),np.cos(DH_alpha[k-1]) * DH_d[k-1]],
			[0,0,0,1]])

		return Tk

	def calc_I(self):

		link_mass = self.link_mass
		link_length = self.link_length
		link_radius = self.link_radius
		num_joints = self.num_joints


		# Inertia tensors arround the centers of mass
		I_com = np.zeros((num_joints,3,3))
		for i in range(0,num_joints):
			tempa = (link_mass[i]*np.power(link_length[i],2)/12) + (link_mass[i]*np.power(link_radius[i],2)/4)
			tempb = link_mass[i]*np.power(link_radius[i],2)/2
			I_com[i] = tempa * np.identity(3)
			if(self.s_hat[i][0] != 0):
				I_com[i][0,0] = tempb
			elif(self.s_hat[i][1] != 0):
				I_com[i][1,1] = tempb
			elif(self.s_hat[i][2] != 0):
				I_com[i][2,2] = tempb

		return I_com

	def calc_I_hat(self):

		num_joints = self.num_joints
		s_hat = self.s_hat
		link_mass = self.link_mass

		I_com = self.calc_I()

		# Inertia tensors around the shifted rotation points
		I_hat = np.zeros((num_joints,3,3))
		s_hat_cross = np.zeros((num_joints,3,3))
		for i in range(0,num_joints):
			s_hat_cross[i] = np.matrix([[0,-s_hat[i][2],s_hat[i][1]],[s_hat[i][2],0,-s_hat[i][0]],[-s_hat[i][1],s_hat[i][0],0]])
			I_hat[i] = I_com[i] - link_mass[i] * np.dot(s_hat_cross[i],s_hat_cross[i])

		return I_hat


	def calc_H_hat(self):

		link_mass = self.link_mass
		s_hat = self.s_hat
		num_joints = self.num_joints

		I_hat = self.calc_I_hat()

		# Calculate the H_hat matrix
		self.H_hat = np.zeros((num_joints,4,4))
		for i in range(0,num_joints):
			self.H_hat[i] = np.matrix([[(-I_hat[i][0][0] + I_hat[i][1][1] + I_hat[i][2][2])/2, 0, 0, link_mass[i]*s_hat[i][0]],
			[0, (I_hat[i][0][0] - I_hat[i][1][1] + I_hat[i][2][2])/2, 0, link_mass[i]*s_hat[i][1]],
			[0, 0, (I_hat[i][0][0] + I_hat[i][1][1] - I_hat[i][2][2])/2, link_mass[i]*s_hat[i][2]],
			[link_mass[i]*s_hat[i][0], link_mass[i]*s_hat[i][1], link_mass[i]*s_hat[i][2], link_mass[i]]])

# Calculate the g matrix
	def calc_g(self,theta):

		s_hat = self.s_hat
		link_mass = self.link_mass

		g_tilda = [0,0,-0.0981,0]
		g = np.zeros(self.num_joints)
		for i in range(0,self.num_joints):
			for j in range(0,self.num_joints):
				T_ji = self.calc_dTi_dqj(j+1,i+1,theta)
				s_hat_resized = np.concatenate([s_hat[j],[1]])
				g[i] -= link_mass[j] * np.linalg.multi_dot([g_tilda,T_ji,np.transpose(s_hat_resized)])
		return g



	def calc_dTi_dqj(self,i,j,theta):

		assert len(theta) == 5, "Error: calc_dTi_dqj - Input parameter theta should be a vector."
		assert (i >= 1) and (i <= 5), "Error: calc_dTi_dqj - indexes starting with one."
		assert (j >= 1) and (j <= 5), "Error: calc_dTi_dqj - indexes starting with one."
		#assert (j <= i), "Error: calc_dTi_dqj - j must be less than or equal to i."
		
		if j > i :
			return np.zeros((4,4))
		T = np.identity(4)
		delta = np.matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
		for n in range(0,i):
			T = T.dot(self.calc_Tk(n+1,theta[n]))
			#print("In here")
			if n == j-1:
				T = T.dot(delta)
				#print("Delta Pow")
		return T

		# Calculate the h matrix
	def calc_h(self,theta,theta_dot):
		h = np.zeros(self.num_joints)
		for i in range(0,self.num_joints):

			for j in range(0,self.num_joints):
				for m in range(0,self.num_joints):
					for k in range(max(i,j,m),self.num_joints):
						T_kjm = self.calc_dTi_dqj_dqk(k+1,j+1,m+1,theta)
						T_ki = self.calc_dTi_dqj(k+1,i+1,theta)

						h[i] += np.trace(np.linalg.multi_dot([T_kjm,self.H_hat[k],np.transpose(T_ki)]))*theta_dot[j]*theta_dot[m]
	

		return h


	#i,j,k --> k,j.m
	def calc_dTi_dqj_dqk(self,i,j,k,theta):

		assert len(theta) == 5, "Error: calc_dTi_dqj - Input parameter theta should be a vector."
		assert (i >= 1) and (i <= 5), "Error: calc_dTi_dqj - indexes starting with one."
		assert (j >= 1) and (j <= 5), "Error: calc_dTi_dqj - indexes starting with one."
		assert (k >= 1) and (k <= 5), "Error: calc_dTi_dqj - indexes starting with one."
		assert (j <= i), "Error: calc_dTi_dqj - j must be less than or equal to i."
		

		# print("===== "," i: ",i," j: ",j," k: ",k)

		T = []

		if (max(j,k) > i):
			T = np.matrix(np.zeros((4,4)))
		elif (i >= k and i >= j):
			delta = np.matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
			T = np.identity(4)
			for n in range(0,i):
				T = T.dot(self.calc_Tk(n+1,theta[n]))
				if n == j-1:
					T = T.dot(delta)
		else:
			assert (3 > 4), "Error: calc_dTi_dqj_dqk - some kind of problem here."


		return T

	def calc_M(self,theta):
		M = np.zeros((self.num_joints,self.num_joints))
		for i in range(0,self.num_joints):
			for j in range(0,self.num_joints):
				for k in range(max(i,j),self.num_joints):
					T1 = self.calc_dTi_dqj(k+1,j+1,theta)
					T2 = np.transpose(self.calc_dTi_dqj(k+1,i+1,theta))
					M[i][j] += np.trace(np.linalg.multi_dot([T1,self.H_hat[k],T2]))
		return M


	def run_dynamics(self,i,torque_in,dt=1):

		M = self.calc_M(self.q[i])
		M_inv = np.linalg.inv(M)

		self.q_ddot[i] =  np.dot(M_inv,torque_in)

		self.q_dot[i+1] = self.q_dot[i] + self.q_ddot[i]*dt
		self.q[i+1] = self.q[i] + self.q_dot[i+1]*dt


	def loop_dynamics_with_g(self,i,torque,dt=1):

		M = self.calc_M(q[pos])
		g = self.calc_g(q[pos])
		M_inv = np.linalg.inv(M)
		q_ddot[pos] =  np.dot(M_inv,(torque-g))
		q_dot[pos+1] = q_dot[pos] + q_ddot[pos] * time_step
		q[pos+1] = q[pos] + q_dot[pos+1] * time_step

		return num_iter,time_step,q,q_dot,q_ddot

	def loop_dynamics_with_g_h_loop(self,q,q_dot,q_ddot,pos,lines,torque,time_step):

		# Loop through iterations
		M = self.calc_M(q[pos])
		h = self.calc_h(q[pos],q_dot[pos])
		g = self.calc_g(q[pos])
		M_inv = np.linalg.inv(M)

		q_ddot[pos] =  np.dot(M_inv,(torque-h-g))
		q_dot[pos+1] = q_dot[pos] + q_ddot[pos] * time_step
		q[pos+1] = q[pos] + q_dot[pos+1] * time_step


	def convert_coordinate_frames(self,i,scale_factor=0.05,endpoint=True): 

		link_coordinates = np.zeros((self.num_joints+1,4))
		link_coordinates[:,3] = 1
		link_coordinates[1][2] = -15

		q_base_reference = np.zeros((self.num_joints+1,4))
		q_base_reference[:,3] = 1


		T0i = self.calc_Tk(1,self.q[i][0])
		q_base_reference[1,:] = np.matmul(T0i,link_coordinates[1,:])

		print(q_base_reference[1,:])

		for j in range(2,self.num_joints+1):


			#calculate transformation matrix for each link frame
			T0i = np.matmul(T0i,self.calc_Tk(j,self.q[i][j-1]))

			#transform the link coordinates to the correct reference frame
			q_base_reference[j,:] = np.dot(T0i,link_coordinates[j,:])

			print(q_base_reference[j,:])

		if(endpoint):
			end_peice = np.matmul(T0i,[10,10,0,1])
			q_base_reference = np.concatenate([q_base_reference,end_peice],axis=0)


		return q_base_reference*scale_factor




	







