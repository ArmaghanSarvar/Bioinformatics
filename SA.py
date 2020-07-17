import numpy as np
import math


class Node():
    def __init__(self, value, parent):
        self.value = value
        self.parent = parent


match_award = 1
mismatch_penalty = -1
open_gap_penalty = -10
extend_gap_penalty = -1


def get_max_M(row, col, S):
    num1 = M[row - 1][col - 1].value + S
    num2 = Ix[row - 1][col - 1].value + S
    num3 = Iy[row - 1][col - 1].value + S
    value = max(num1, num2, num3)
    parent = []
    if value == num1:
        parent.append([row - 1, col - 1, 'M'])
    if value == num2:
        parent.append([row - 1, col - 1, 'X'])
    if value == num3:
        parent.append([row - 1, col - 1, 'Y'])
    return Node(value, parent)


def get_max_Ix(row, col):
    num1 = M[row - 1][col].value + open_gap_penalty
    num2 = Ix[row - 1][col].value + extend_gap_penalty
    value = max(num1, num2)
    parent = []
    if value == num1:
        parent.append([row - 1, col, 'M'])
    if value == num2:
        parent.append([row - 1, col, 'X'])
    return Node(value, parent)


def get_max_Iy(row, col):
    num1 = M[row][col - 1].value + open_gap_penalty
    num2 = Iy[row][col - 1].value + extend_gap_penalty
    value = max(num1, num2)
    parent = []
    if value == num1:
        parent.append([row, col - 1, 'M'])
    if value == num2:
        parent.append([row, col - 1, 'Y'])
    return Node(value, parent)


s1 = input("Enter The First Sequence : ")
s2 = input("Enter The Second Sequence : ")

rows = len(s1) + 1
cols = len(s2) + 1

# init M
M = np.full((rows, cols), Node)
M[0][0] = Node(0, [])

for i in range(1, rows):
    M[i][0] = Node(-math.inf, [[i - 1, 0, 'M']])
for i in range(1, cols):
    M[0][i] = Node(-math.inf, [[0, i - 1, 'M']])

# init Ix
Ix = np.full((rows, cols), Node)
Ix[0][0] = Node(open_gap_penalty + (0 - 1) * extend_gap_penalty, [])

for i in range(1, rows):
    Ix[i][0] = Node(open_gap_penalty + (i - 1) * extend_gap_penalty, [[i - 1, 0, 'X']])

for i in range(1, cols):
    Ix[0][i] = Node(-math.inf, [[0, i - 1, 'X']])


# init Iy
Iy = np.full((rows, cols), Node)
Iy[0][0] = Node(open_gap_penalty + (0 - 1) * extend_gap_penalty, [])

for i in range(1, rows):
    Iy[i][0] = Node(-math.inf, [[i - 1, 0, 'Y']])

for i in range(1, cols):
    Iy[0][i] = Node(open_gap_penalty + (i - 1) * extend_gap_penalty, [[0, i - 1, 'Y']])

for row in range(1, rows):
    for col in range(1, cols):
        # match
        if s1[row - 1] == s2[col - 1]:
            S = match_award
            M[row][col] = get_max_M(row, col, S)
        else:
            S = mismatch_penalty
            M[row][col] = get_max_M(row, col, S)
        Ix[row][col] = get_max_Ix(row, col)
        Iy[row][col] = get_max_Iy(row, col)


last1 = M[rows - 1][cols - 1].value
last2 = Ix[rows - 1][cols - 1].value
last3 = Iy[rows - 1][cols - 1].value
start = max(last1, last2, last3)
start_nodes = []
start_node = Node

if start == last1:
    start_node = M[rows - 1][cols - 1]
    start_nodes.append(start_node)

if start == last2:
    start_node = Ix[rows - 1][cols - 1]
    start_nodes.append(start_node)

if start == last3:
    start_node = Iy[rows - 1][cols - 1]
    start_nodes.append(start_node)



cur_row = rows - 1
cur_col = cols - 1

score = start_node.value
def get_gapped(res):
    min_gap = len(res) - 1
    min_gap_idx = len(res[0]) - 1
    i = 0
    while i < len(res):
        for j in range(len(res[0])):
            if res[i][j] == '-' or res[i + 1][j] == '-':
                if j <= min_gap_idx:
                    min_gap_idx = j
                    min_gap = i
        i += 2
    res.remove(res[min_gap])
    res.remove(res[min_gap])
    return res


def traceback(current_node, c1, c2, current_col, current_row):
    if not current_node.parent:
        return [c1, c2]
    paths = []
    for parent in current_node.parent:
        col = parent[1]
        row = parent[0]
        if parent[2] == 'M':
            neighbor = M[parent[0], parent[1]]
        elif parent[2] == 'X':
            neighbor = Ix[parent[0], parent[1]]
        elif parent[2] == 'Y':
            neighbor = Iy[parent[0], parent[1]]
        if row == current_row:
            paths += traceback(neighbor, c1 + '-', c2 + s2[col], col, row)
        elif col == current_col:
            paths += traceback(neighbor, c1 + s1[row], c2 + '-', col, row)
        else:
            paths += traceback(neighbor, c1 + s1[row], c2 + s2[col], col, row)
    return paths


final = []
for start_node in start_nodes:
    results = traceback(start_node, "", "", cur_col, cur_row)
    for i in range(len(results)):
        results[i] = results[i][::-1]
    final.append(results)

final_final = []
for i in range(len(final)):
    for j in range(len(final[i])):
        final_final.append(final[i][j])


while True:
    if len(final_final) == 2:
        break
    final_final = get_gapped(final_final)


print('Final Alignment is: ', final_final)
print('Final Score is: ', score)
