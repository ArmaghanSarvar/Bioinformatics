import numpy as np
import math


class Node:
    def __init__(self, value, parent):
        self.value = value
        self.parent = parent


match_award = 1
mismatch_penalty = -1
open_gap_penalty = -10
extend_gap_penalty = -1


def get_max_M(row, col, S, M, Ix, Iy):
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


def get_max_Ix(row, col, M, Ix):
    num1 = M[row - 1][col].value + open_gap_penalty
    num2 = Ix[row - 1][col].value + extend_gap_penalty
    value = max(num1, num2)
    parent = []
    if value == num1:
        parent.append([row - 1, col, 'M'])
    if value == num2:
        parent.append([row - 1, col, 'X'])
    return Node(value, parent)


def get_max_Iy(row, col, M, Iy):
    num1 = M[row][col - 1].value + open_gap_penalty
    num2 = Iy[row][col - 1].value + extend_gap_penalty
    value = max(num1, num2)
    parent = []
    if value == num1:
        parent.append([row, col - 1, 'M'])
    if value == num2:
        parent.append([row, col - 1, 'Y'])
    return Node(value, parent)


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


def traceback(current_node, c1, c2, current_col, current_row, M, Ix, Iy, s1, s2):
    if not current_node.parent:
        return [c1, c2]
    paths = []
    for parent in current_node.parent:
        row = parent[0]
        col = parent[1]
        if parent[2] == 'M':
            neighbor = M[parent[0], parent[1]]
        elif parent[2] == 'X':
            neighbor = Ix[parent[0], parent[1]]
        elif parent[2] == 'Y':
            neighbor = Iy[parent[0], parent[1]]
        if row == current_row:
            paths += traceback(neighbor, c1 + '-', c2 + s2[col], col, row, M, Ix, Iy, s1, s2)
        elif col == current_col:
            paths += traceback(neighbor, c1 + s1[row], c2 + '-', col, row, M, Ix, Iy, s1, s2)
        else:
            paths += traceback(neighbor, c1 + s1[row], c2 + s2[col], col, row, M, Ix, Iy, s1, s2)
    return paths


def SA(s1, s2):
    rows = len(s1) + 1
    cols = len(s2) + 1

    M = np.full((rows, cols), Node)
    M[0][0] = Node(0, [])

    for i in range(1, rows):
        M[i][0] = Node(-math.inf, [[i - 1, 0, 'M']])
    for i in range(1, cols):
        M[0][i] = Node(-math.inf, [[0, i - 1, 'M']])

    Ix = np.full((rows, cols), Node)
    Ix[0][0] = Node(open_gap_penalty + (0 - 1) * extend_gap_penalty, [])

    for i in range(1, rows):
        Ix[i][0] = Node(open_gap_penalty + (i - 1) * extend_gap_penalty, [[i - 1, 0, 'X']])

    for i in range(1, cols):
        Ix[0][i] = Node(-math.inf, [[0, i - 1, 'X']])

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
                M[row][col] = get_max_M(row, col, S, M, Ix, Iy)
            else:
                S = mismatch_penalty
                M[row][col] = get_max_M(row, col, S, M, Ix, Iy)
            Ix[row][col] = get_max_Ix(row, col, M, Ix)
            Iy[row][col] = get_max_Iy(row, col, M, Iy)

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

    final = []
    for start_node in start_nodes:
        results = traceback(start_node, "", "", cur_col, cur_row, M, Ix, Iy, s1, s2)
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
    return [final_final, score]


# n = int(input("Enter the Number of Sequences : "))
n = 4
sequences_matrix = np.zeros((n, n))
sequences_list = []
for i in range(n):
    # s = input("Enter The " + str(i + 1) + "th Sequence : ")
    s = input()
    sequences_list.append(s)

for i in range(n):
    for j in range(i+1, n):
        score = SA(sequences_list[i], sequences_list[j])[1]
        sequences_matrix[i][j] = score
        sequences_matrix[j][i] = score

# print('The Score Matrix is: ')
# print(sequences_matrix)
row_sum = np.sum(sequences_matrix, axis=1)
center_idx = np.where(row_sum == max(row_sum))[0][0]
center = sequences_list[center_idx]
# print('The Center Sequence is: ', center)
alignments = []
for seq_idx in range(len(sequences_list)):
    if seq_idx != center_idx:
        # print('center is', center)
        # print('seq is', seq)
        res = SA(center, sequences_list[seq_idx])
        alignments.append(res[0])
        # print('Alignment is: ', res[0][0], ' *** ', res[0][1])
        # print('Alignment Score is: ', res[1])


def gap_all_in_MSA(multi, pos):
    for idx in range(1, len(multi)):
        multi[idx] = multi[idx][:pos] + '-' + multi[idx][pos:]
    return multi


def insert_gap_in_MSA(multi):
    for idx in range(1, len(multi)):
        multi[idx] += '-'
    return multi


# print('Alignments are: ', alignments)


MSA = alignments[0]
for i in range(1, len(alignments)):
    center1 = MSA[0]
    # print("center1:", center1)
    center2 = alignments[i][0]
    # print("center2:", center2)
    k1 = 0
    k2 = 0
    while True:
        if k1 == len(center1) - 1:
            while k2 < len(center2) - 1:
                center1 += '-'
                MSA = insert_gap_in_MSA(MSA)
                MSA[0] = center1
                k2 += 1
            break
        if k2 == len(center2) - 1:
            while k1 < len(center1) - 1:
                center2 += '-'
                alignments[i][1] += '-'
                k1 += 1
            break
        if center1[k1] == center2[k2]:
            pass
        if center1[k1] == '-' and center2[k2] != '-':
            center2 = center2[:k2] + '-' + center2[k2:]
            alignments[i][1] = alignments[i][1][:k2] + '-' + alignments[i][1][k2:]
        elif center2[k2] == '-' and center1[k1] != '-':
            center1 = center1[:k1] + '-' + center1[k1:]
            MSA = gap_all_in_MSA(MSA, k1)
            MSA[0] = center1
        k1 += 1
        k2 += 1
    # print("center1:", center1)
    # print("center2:", center2)
    # print("alignments:", alignments)
    MSA.append(alignments[i][1])
    # print("MSA:", MSA)


def omit_gaps(string):
    i = 0
    while i < len(string):
        if string[i] == '-':
            string = string[:i] + string[i + 1:]
            i -= 1
        i += 1
    return string


def find_in_input(s):
    for i in range(len(sequences_list)):
        if sequences_list[i] == s:
            return i


output = {}
for s in MSA:
    s1 = omit_gaps(s)
    output[s] = find_in_input(s1)

sorted = sorted(output.items(), key=lambda x: x[1])

for key, value in sorted:
    print(key)
