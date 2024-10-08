import sympy as sp
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


def sixteen_cell():
    # Diagram:
    #     o2
    # 0  1|   3
    # o---o---o

    # Matrix:
    matrix = sp.Matrix([
        [1, 3, 2, 2],
        [3, 1, 3, 3],
        [2, 3, 1, 2],
        [2, 3, 2, 1],
    ])
    return matrix


def cube():
    matrix = sp.Matrix([
        [1, 4, 2],
        [4, 1, 3],
        [2, 3, 1],
    ])
    return matrix


def triangle():
    matrix = sp.Matrix([
        [1, 3],
        [3, 1],
    ])
    return matrix


def five_cube():
    matrix = sp.Matrix([
        [1, 4, 2, 2, 2],
        [4, 1, 3, 2, 2],
        [2, 3, 1, 3, 2],
        [2, 2, 3, 1, 3],
        [2, 2, 2, 3, 1],
    ])
    return matrix


def reflect(d, n):
    return d - 2 * d.dot(n) * n


def reflection_matrix(n):
    eye = sp.eye(N)
    return sp.Matrix.hstack(
        *[reflect(eye.col(i), n) for i in range(N)]
    )


# ==== Initial Coxeter matrix ====
matrix = sixteen_cell()
N = sp.shape(matrix)[0]

print(N, matrix)


# ==== Generate normal vectors ====
normals = sp.zeros(1, N)
normals[0] = 1


for i in range(1, N):
    subarray = normals[:i, :i]

    b = matrix[:i, i]
    b = b.applyfunc(lambda x: sp.cos(sp.pi / x))

    solution = subarray.inv() * b

    new_normal = sp.zeros(1, N)

    for j in range(i):
        new_normal[0, j] = solution[j]

    sum_of_squares = new_normal.dot(new_normal)

    last_element = sp.sqrt((1 - sum_of_squares))

    new_normal[i] = last_element

    row_idx = normals.shape[0]
    normals = normals.row_insert(row_idx, new_normal)


# ==== Ocular validation of normal vectors ====
for i in range(N):
    for j in range(i, N):
        norm_i = normals[i, :]
        norm_j = normals[j, :]

        dp = norm_i.dot(norm_j)
        if dp == 1:
            coxeter_number = 1
        else:
            coxeter_number = 1 / (sp.acos(dp) / sp.pi)
        print((i, j), coxeter_number, matrix[i, j])
print("normals:", normals)

# ==== Generators of Coxeter group (as matrices) ====
gens = []
for i in range(N):
    n = normals.row(i).T
    gens.append(reflection_matrix(n))

for g in gens:
    print(g)


class GroupElement:
    sequence: list[int]
    matrix: sp.ImmutableMatrix

    def __init__(self, sequence, matrix):
        self.sequence = sequence
        self.matrix = matrix

    def __eq__(self, o):
        return self.matrix == o.matrix

    def __hash__(self):
        return hash(self.matrix)


# ==== Get all elements of Coxeter group ====
elements = []

queue = [
    GroupElement([], sp.ImmutableMatrix(sp.eye(N)))
]
seen = set(queue)

while queue:
    ele = queue.pop(0)

    elements.append(ele)

    for i, g in enumerate(gens):
        sequence = ele.sequence + [i]
        matrix = g * ele.matrix
        new_ele = GroupElement(sequence, matrix)

        if new_ele not in seen:
            queue.append(new_ele)
            seen.add(new_ele)

# Check order
print(len(elements))

root_point = sp.Matrix([2, 2, 1])

all_points = []
for ele in elements:
    point = ele.matrix * root_point
    all_points.append(point.T)

all_points = sp.Matrix(all_points)

hull = ConvexHull(all_points)
print(hull.simplices)
print(hull.simplices.shape)


# plottttt
ax = plt.figure().add_subplot(projection='3d')
ax.scatter(all_points[:, 0], all_points[:, 1], all_points[:, 2])
plt.show()
