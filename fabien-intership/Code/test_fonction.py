import numpy as np
import matplotlib.pyplot as plt

matrice = np.array([[1,1,2],[2,23,6],[8,9,15]])
weight = np.array([[2],[3],[4]])
test = np.array([[2]])
print(matrice[:, 0])
print(matrice[:, 0].reshape(3,1))

matrice[:, 0] = weight.flatten()
print(matrice)

calcul = matrice @ weight
print(calcul)
print(weight[0,0])
print(test.T)