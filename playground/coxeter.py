import numpy as np

# Diagram:
#     o2
# 0  1|   3
# o---o---o

# Matrix:
matrix = np.array([
    [1, 3, 2, 2],
    [3, 1, 3, 3],
    [2, 3, 1, 2],
    [2, 3, 2, 1],
])
N = len(matrix)

normals = np.array([[1, 0, 0, 0]])


for i in range(1, N):
    subarray = normals[:i, :i]
    b = np.cos(np.pi / matrix[i, :i])
    print(subarray)
    print(b)

    solution = np.linalg.inv(subarray) @ b
    print(solution)

    new_normal = np.zeros(N)
    new_normal[:i] = solution

    print(new_normal)
    last_element = (1 - np.sum(new_normal**2)) ** 0.5
    new_normal[i] = last_element
    print(new_normal)
    print(np.sum(new_normal**2))

    normals = np.vstack((normals, [new_normal]))

    print()


for i in range(N):
    for j in range(i, N):
        norm_i = normals[i]
        norm_j = normals[j]

        dp = norm_i @ norm_j
        if dp == 1:
            coxeter_number = 1
        else:
            coxeter_number = 1 / (np.acos(dp) / np.pi)
        print((i, j), coxeter_number, matrix[i, j])
print(normals)
