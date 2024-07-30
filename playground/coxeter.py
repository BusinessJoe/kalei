import sympy as sp


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


def reflect(d, n):
    return d - 2 * d.dot(n) * n


def reflection_matrix(n):
    x = reflect(sp.Matrix([1, 0, 0]), n)
    y = reflect(sp.Matrix([0, 1, 0]), n)
    z = reflect(sp.Matrix([0, 0, 1]), n)
    return x.row_join(y).row_join(z)


matrix = cube()
N = sp.shape(matrix)[0]

print(N, matrix)

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

gens = []
for i in range(N):
    n = normals.row(i).T
    gens.append(reflection_matrix(n))

for g in gens:
    print(g)


# class GroupElement:
#     sequence: list[int]
#     matrix: sp.Matrix

#     def __init__(self, sequence: list[int], matrix: sp.Matrix):
#         self.sequence = sequence
#         self.matrix = matrix


# elements = []
# queue = [
#     GroupElement([], sp.eye(N))
# ]

# while queue:

