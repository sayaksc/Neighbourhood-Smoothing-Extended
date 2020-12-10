import numpy as np
import math
import random
import matplotlib.pyplot as plt
from graphons import graphon1 as graphon


def nbdsmooth(A, A1, A2, D):
	n = len(A)
	h = math.sqrt(math.log(n)*1.0/n)
	kernel_mat = np.zeros((n, n))
	for i in range(n):
		kernel_mat[i] =  (D[i]<np.quantile(D[i], h)).astype(float)
	kernel_mat = kernel_mat/(np.reshape(np.repeat(np.sum(kernel_mat, 1), n), (n, n)) + 0.1)
	P = kernel_mat.dot(A)
	P = (P + np.transpose(P))/2
	return P





def D_ext2(n, A1, A2):
	n1 = len(A1)
	n2 = len(A2)

	A1sq = np.dot(A1, A1)
	A1sq = A1sq - np.diag(np.diag(A1sq))
	A2sq = np.dot(A2, A2)
	
	A2sq = A2sq - np.diag(np.diag(A2sq))
	D1 = np.zeros((n1, n1))
	D2 = np.zeros((n2, n2))
	Dtemp1 = np.zeros((n, n))
	Dtemp2 = np.zeros((n, n))
	D = np.zeros((n, n))

	for i in range(n1):
		for j in range(n1):
			vi = np.array(A1sq[i])
			vj = np.array(A1sq[j])
			vi[j] = 0
			vj[i] = 0
			D1[i, j] = np.max(abs(vi - vj))*1.0 / n1


	for i in range(n2):
		for j in range(n2):
			vi = np.array(A2sq[i])
			vj = np.array(A2sq[j])
			vi[j] = 0
			vj[i] = 0
			D2[i, j] = np.max(abs(vi - vj))*1.0 / n2

	Dtemp1[:n1, :n1] = np.array(D1)
	Dtemp2[-n2:, -n2:] = np.array(D2)
	D = Dtemp1 + Dtemp2
	D[(n - n2):n1, (n - n2):n1] = 0.5*D[(n - n2):n1, (n - n2):n1]
	for i in range(n-n2):
		for j in range(n1, n):
			temp1 = D[i, (n - n2):n1] + D[j, (n - n2):n1]
			temp2 = abs(D[i, (n - n2):n1] - D[j, (n - n2):n1])
			temp3 = np.min(temp1)
			temp4 = np.max(temp2)
			# D[i, j] = (temp3 + temp4)/2.0
			# D[j, i] = (temp3 + temp4)/2.0
			# temp5 = (temp3 + temp4)/2.0
			temp5 = 2.0*temp3*temp4/(temp3 + temp4)
			# temp5 = math.sqrt(temp3*temp4*1.0)
			# temp5 = math.sqrt((temp3*temp3 + temp4*temp4)/2.0)
			# temp5 = temp3
			# temp5 = temp4
			D[i, j] = temp5
			D[j, i] = temp5
	return D

def read_file(f, n):
        file = open(f, "r")
        N = np.zeros((n, n))
        for i in range(n):
                L = file.readline().split()
                for j in range(n):
                        N[i][j] = float(L[j])
        return N



def write_file(f, P, n):
        file = open(f, "w")
        for i in range(n):
                s = ""
                for j in range(n):
                        s += str(float(P[i][j])) + " "
                s += "\n"
                file.write(s)
        file.close()
        return
def graphon(x, y, K = 2, p = 0.3, q = 0.03):
    # return (math.sin(math.pi*5*(x+ y -1) + 1)/2+0.5)
    # return 1-max(x, y)
    # return (1 - (1 + math.exp(-15*(0.8*abs(x-y))**(4/5) - 0.1))**(-1))
    # return ((x**2 + y**2)/3)*math.cos(1/(x**2 + y**2)) + 0.15
    # return 1/(1 + math.exp(-x -y))
    t1 = int(K*x)
    t2 = int(K*y)
    if (t1 == t2):
        return p
    else:
        return q


n = 1000
L = np.random.uniform(0, 1, n)
perm = np.argsort(L)
P_est = np.zeros((n, n))
A = np.zeros((n, n))
for i in range(n):
	for j in range(i+1, n):
		x = L[i]
		y = L[j]
		P_est[i, j] = graphon(x, y)
		P_est[j, i] = P_est[i, j]
		p = random.random()
		if p < P_est[i, j]:
			A[i, j] = 1
			A[j, i] = 1

error_fro = []
error_max = []
total = []
# M = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.499]
M = [0.05, 0.15, 0.3]
# M = [0.49]
for m in M:
	L_fro = []
	L_max = []
	print("--------------", m, "------------------")
	n1 = int((0.5+m)*n)
	n2 = int((0.5+m)*n)
	B = np.zeros((n, n))

	#Create Matrix
	for i in range(n):
		for j in range(i, n):
			if (i<n1 and j<n1) or (n-i<n2-1 and n-j<n2-1) :
				B[i][j] = A[i][j]
				B[j][i] = A[i][j]


	#Check the partially revealed graph
	# if (m == 0.1):
	# 	temp = 1-np.array(B)
	# 	temp[0, 0] = -1
	# 	temp[:(n-n1), n1:] = -1
	# 	temp[n1:, :(n-n2)] = -1
	# 	temp2 = np.zeros((n, n))
		# for i in range(n):
		# 	for j in range(n):
		# 		temp2[i, j] = temp[perm[i], perm[j]]
		# plt.imshow(temp, cmap='hot', interpolation='nearest', origin = "lower")
		# plt.title("Partially Revealed")
		# plt.savefig("../datasets/" + string + "/covered.png")
		# plt.clf()
		
		# plt.show()

	#Calculate distance matrix
	E = D_ext2(n, B[:n1, :n1], B[-n2:, -n2:])

	#perform NBS by finding neighbours everytime
	P = nbdsmooth(B, B[:n1, :n1], B[-n2:, -n2:], E)
	P1 = np.zeros((n, n))
	p_fro = np.linalg.norm((P1-P), ord = "fro")/(n*1.0)
	p_max = np.max(np.abs(P1-P))
	P1 = np.zeros((n, n))


	while(abs(p_fro) > 0.001):

	# for itern in range(10):
		print(p_fro, p_max)
		B = np.array(P)
		P1 = np.array(P)
		for i in range(n):
			for j in range(i, n):
				if (i<n1 and j<n1) or (n-i<n2-1 and n-j<n2-1) :
					B[i][j] = A[i][j]
					B[j][i] = A[i][j]


		# D = D_ext2(n, B[:n1, :n1], B[-n2:, -n2:])
		P = nbdsmooth(B, B[:n1, :n1], B[-n2:, -n2:], E)
		# p_err = np.linalg.norm((P_est-P), ord = "fro")/(n*1.0)
		p_fro = np.linalg.norm((P1-P), ord = "fro")/(n*1.0)
		p_max = np.max(np.abs(P1-P))
		L_max.append(p_max)
		L_fro.append(p_fro)
		# p = np.linalg.norm((P1-P), ord = "fro")/(n*1.0)
	# print("list:", L_rr)
	# P = nbdsmooth(B, B[:n1, :n1], B[-n2:, -n2:], E)
	error_fro.append(L_fro)
	error_max.append(L_max)


	# T = np.zeros((n, n))

	# for i in range(n):
		# for j in range(n):
			# T[i, j] = P[perm[i], perm[j]]
	# plt.imshow(1-T, cmap='gray', interpolation='nearest', origin = "lower")
	# plt.title("$m = " + str(m) + "n$")
	# plt.show()

	


# T = np.zeros((n, n))

# for i in range(n):
# 	for j in range(n):
# 		T[i, j] = A[perm[i], perm[j]]


# T[0, 0] = -1
# plt.imshow(T, cmap='hot', interpolation='nearest', origin = "lower")
# plt.title("Graph")
# plt.savefig("../datasets/" + string + "/graph.png")
# plt.clf()
# plt.show()


#inverse permutation
# T = np.zeros((n, n))

# for i in range(n):
# 	for j in range(n):
# 		T[i, j] = P_est[perm[i], perm[j]]
# plt.imshow(1-T, cmap='gray', interpolation='nearest', origin = "lower")
# plt.title("Graphon")
# plt.show()


c10 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.plot(c10, error_fro[0], label = M[0], color = "red")
plt.plot(c10, error_max[0], "--", color = "red")

plt.plot(c10, error_fro[1], label = M[1], color = "green")
plt.plot(c10, error_max[1], "--",  color = "green")

plt.plot(c10, error_fro[2], label = M[2], color = "blue")
plt.plot(c10, error_max[2], "--",  color = "blue")

plt.legend()
plt.show()