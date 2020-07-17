import os
import numpy as np
import time
import sys

DNA_list = ['A', 'T', 'C', 'G', '-']
amino_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']

pseudo_count = 0.001


def is_DNA(sequence):
    for s in range(len(sequence)):
        if sequence[s] != 'A' and sequence[s] != 'T' and sequence[s] != 'C' and sequence[s] != 'G' and sequence[s] != '-':
            return False
    return True


def read_dir():
    if len(sys.argv) < 2:
        path = input('Enter dir of test cases: ')
    else:
        path = sys.argv[1]
    onlyfiles = next(os.walk(path + '/in/'))[2]
    return len(onlyfiles), path


def initialization(path, file_number):
    global insert_states, sequences_bef_update
    input_path = path + '/in/'
    output_path = path + '/out/'

    with open(input_path + 'input' + str(file_number + 1) + '.txt') as f:
        my_list = list(f)
    n = int(my_list[0].split(" ")[0])
    gap_threshold = int(my_list[0].split(" ")[1])
    test_cases = []
    length = 0
    sequences_list = []
    for i in range(n):
        s = my_list[i + 1]
        length = len(s) - 1
        sequences_list.append(s[:-1])
    test_cases.append(my_list[-1][:-1])
    with open(output_path + 'output' + str(file_number + 1) + '.txt') as f:
        my_list2 = list(f)
        true_output = my_list2[0][:-1]
    sequences_bef_update = sequences_list.copy()
    gaps = [0] * length
    insert_states = [0] * length
    for i in range(n):
        for j in range(length):
            if sequences_list[i][j] == '-':
                gaps[j] += 1

    for i in range(len(gaps)):
        if gaps[i] > gap_threshold:
            insert_states[i] = 1

    i = 0
    boolean = False
    while i < len(gaps):
        if gaps[i] > gap_threshold:
            boolean = True
            for j in range(n):
                sequences_list[j] = sequences_list[j][:i] + sequences_list[j][i + 1:]
            i -= 1
        i += 1
        if boolean:
            gaps.remove(gaps[i])
            boolean = False
    # print("Sequences:\n", sequences_list)
    # print(insert_states)
    return sequences_list, is_DNA(sequences_list[0]), test_cases, true_output


def init_emission(emission_matrix, is_DNA, sequences_list):
    ns = len(sequences_list)
    residue_list = DNA_list if is_DNA else amino_list
    delete_counter = 0
    del_res_counter = 0
    match_counter = 0
    match_res_counter = 0

    for state_idx in range(emission_matrix.shape[0]):
        for residue_idx in range(len(residue_list)):
            if state_idx == 0 or state_idx == emission_matrix.shape[0] - 1:
                emission_matrix[state_idx][residue_idx] = 0

            elif state_idx % 3 == 0:  # Deletion
                if del_res_counter % len(residue_list) == 0:
                    column = np.array([sequences_list[j][delete_counter] for j in range(ns)])
                    delete_counter += 1
                    number_of_gaps = np.where(column == "-")[0].shape[0]
                    emission_matrix[state_idx][emission_matrix.shape[1] - 1] = 1
                del_res_counter += 1

            elif state_idx % 3 == 1:  # Insert
                if residue_idx != len(residue_list) - 1:
                    emission_matrix[state_idx][residue_idx] = 1/4 if is_DNA else 1/20

            elif state_idx % 3 == 2:  # Match
                column = np.array([sequences_list[j][match_counter] for j in range(ns)])
                emission_prob = 0
                for r in column:
                    if r == residue_list[residue_idx]:
                        emission_prob += 1
                if residue_idx != emission_matrix.shape[1] - 1:
                    if emission_prob != 0:
                        emission_matrix[state_idx][residue_idx] = emission_prob
                    else:
                        emission_matrix[state_idx][residue_idx] = pseudo_count
                    match_res_counter += 1
                if match_res_counter % len(residue_list) == 0:
                    match_counter += 1
    # Normalize Match States
    for state_idx in range(emission_matrix.shape[0]):
        summation = 0
        if state_idx % 3 == 2:  # Match
            for residue_idx in range(len(residue_list)):
                    summation += emission_matrix[state_idx][residue_idx]
            for residue_idx in range(len(residue_list)):
                if summation != 0:
                    emission_matrix[state_idx][residue_idx] = emission_matrix[state_idx][residue_idx]/summation

    # printing(emission_matrix)
    return emission_matrix


def find_column(profstate2, num_of_states):
    index = int(profstate2[1:])
    state = profstate2[0]
    delete_counter = 0
    insert_counter = 0
    match_counter = 0
    for i in range(num_of_states):
        if state == 'D':
            if i % 3 == 0 and i != 0:
                if delete_counter == index:
                    return i
                delete_counter += 1
        elif state == 'I':
            if i % 3 == 1:
                if insert_counter == index:
                    return i
                insert_counter += 1
        else:
            if i % 3 == 2:
                if match_counter == index:
                    return i
                match_counter += 1


def printing(matrix):
    for i in range(matrix.shape[0]):
        print(i, ':    ', end="")
        for j in range(matrix.shape[1]):
            print(matrix[i][j], end =",  ")
        print('\n')


def init_transition(profHMM, transition_matrix):
    delete_counter = 0
    del_res_counter = 1
    insert_counter = 0
    ins_res_counter = 1
    match_counter = 0
    match_res_counter = 1
    num_of_states = transition_matrix.shape[0]
    for state1_idx in range(num_of_states):
        for state2_idx in range(num_of_states):
            # last row and first column
            if state1_idx == num_of_states - 1 or state2_idx == 0:
                transition_matrix[state1_idx][state2_idx] = 0
            # start state / first row
            if state1_idx == 0:
                for profstate2 in profHMM['S']:
                    state2_found = find_column(profstate2, num_of_states)
                    transition_matrix[state1_idx][state2_found] = profHMM['S'][profstate2]

            if state1_idx % 3 == 0 and state1_idx != 0:  # Deletion
                for profstate2 in profHMM['D' + str(delete_counter)]:
                    if profstate2 == 'E':
                        transition_matrix[num_of_states - 3][num_of_states - 1] = profHMM['D' + str(delete_counter)][profstate2]
                    else:
                        state2_found = find_column(profstate2, num_of_states)
                        transition_matrix[state1_idx][state2_found] = profHMM['D' + str(delete_counter)][profstate2]
                if del_res_counter % num_of_states == 0:
                    delete_counter += 1
                del_res_counter += 1

            elif state1_idx % 3 == 1:  # Insert
                for profstate2 in profHMM['I' + str(insert_counter)]:
                    if profstate2 == 'E':
                        transition_matrix[num_of_states - 2][num_of_states - 1] = profHMM['I' + str(insert_counter)][profstate2]
                    else:
                        state2_found = find_column(profstate2, num_of_states)
                        transition_matrix[state1_idx][state2_found] = profHMM['I' + str(insert_counter)][profstate2]
                if ins_res_counter % num_of_states == 0:
                    insert_counter += 1
                ins_res_counter += 1

            elif state1_idx % 3 == 2 and state1_idx != (num_of_states - 1):  # Match
                # print('state: ', state1_idx)
                # print('mc ', match_counter)
                for profstate2 in profHMM['M' + str(match_counter)]:
                    if profstate2 == 'E':
                        transition_matrix[num_of_states - 4][num_of_states - 1] = profHMM['M' + str(match_counter)][profstate2]
                    else:
                        state2_found = find_column(profstate2, num_of_states)
                        transition_matrix[state1_idx][state2_found] = profHMM['M' + str(match_counter)][profstate2]
                if match_res_counter % num_of_states == 0:
                    match_counter += 1
                match_res_counter += 1
    # printing(transition_matrix)
    return transition_matrix


def init_pi(num_of_states):
    initial_distribution = [1]
    for i in range(num_of_states - 1):
        initial_distribution.append(0)
    return np.array(initial_distribution)


def convert_to_log(matrix):
    for idx1 in range(matrix.shape[0]):
        for idx2 in range(matrix.shape[1]):
            if matrix[idx1][idx2] != 0:
                matrix[idx1][idx2] = np.abs(np.log2(matrix[idx1][idx2]))
    for state_idx in range(matrix.shape[0]):
        summation = 0
        for idx2 in range(matrix.shape[1]):
            summation += matrix[state_idx][idx2]
        for idx2 in range(matrix.shape[1]):
            if summation != 0:
                matrix[state_idx][idx2] = matrix[state_idx][idx2] / summation
    return matrix


def HMMStructure(sequences_list, is_DNA):
    match_states_num = len(sequences_list[0])
    del_states_num = len(sequences_list[0])
    insert_states_num = len(sequences_list[0]) + 1
    all_states = match_states_num + del_states_num + insert_states_num + 2  # start and end
    transition_matrix = np.zeros((all_states, all_states))
    emission_matrix = np.zeros((all_states, 5)) if is_DNA else np.zeros((all_states, 21))
    emission_matrix = init_emission(emission_matrix, is_DNA, sequences_list)
    profHMM = {}  # Holds the transition probs
    match_ids = []
    insert_ids = []
    delete_ids = []
    for i in range(0, match_states_num):
        profHMM.update({"M" + str(i): {}})
        match_ids.append("M" + str(i))
    for i in range(0, del_states_num):
        profHMM.update({"D" + str(i): {}})
        delete_ids.append("D" + str(i))
    for i in range(0, insert_states_num):
        profHMM.update({"I" + str(i): {}})
        insert_ids.append("I" + str(i))
    profHMM.update({"S": {}})
    profHMM.update({"E": {}})

    # Assign transition possibilities:
    # Insert
    for i in range(0, insert_states_num - 1):
        profHMM[insert_ids[i]].update({insert_ids[i]: 0, match_ids[i]: 0, delete_ids[i]: 0})
    # Delete
    for i in range(0, del_states_num - 1):
        profHMM[delete_ids[i]].update({insert_ids[i + 1]: 0, match_ids[i + 1]: 0, delete_ids[i + 1]: 0})
    # Match
    for i in range(0, match_states_num - 1):
        profHMM[match_ids[i]].update({insert_ids[i + 1]: 0, match_ids[i + 1]: 0, delete_ids[i + 1]: 0})

    profHMM["S"].update({insert_ids[0]: 0, match_ids[0]: 0, delete_ids[0]: 0})

    profHMM[insert_ids[insert_states_num - 1]].update({insert_ids[insert_states_num - 1]: pseudo_count, 'E': pseudo_count})
    profHMM[delete_ids[del_states_num - 1]].update({insert_ids[insert_states_num - 1]: pseudo_count,  'E': pseudo_count})
    profHMM[match_ids[match_states_num - 1]].update({insert_ids[insert_states_num - 1]: pseudo_count, 'E': pseudo_count})

    for sequence in sequences_bef_update:
        path = []
        count = 0
        path.append("S")
        for i in range(len(sequence)):
            if sequence[i] != '-' and insert_states[i] != 1:
                path.append("M" + str(count))
                count += 1
            elif sequence[i] == '-' and insert_states[i] != 1:
                path.append("D" + str(count))
                count += 1
            elif sequence[i] != '-' and insert_states[i] == 1:
                path.append("I" + str(count))
            elif sequence[i] == '-' and insert_states[i] == 1:
                nothing = 1
        path.append("E")
        for i in range(len(path) - 1):
            profHMM[path[i]][path[i + 1]] += 1
        # print('Path is: ', path)

    for state1 in profHMM:
        for state2 in profHMM[state1]:
            if profHMM[state1][state2] == 0:
                profHMM[state1][state2] = pseudo_count

    # convert to probability
    for state1 in profHMM:
        summation = 0
        for state2 in profHMM[state1]:
            summation += profHMM[state1][state2]
        if summation != 0:
            for state2 in profHMM[state1]:
                profHMM[state1][state2] = profHMM[state1][state2] / summation
        else:
            for state2 in profHMM[state1]:
                # print(state1, state2)
                profHMM[state1][state2] = 0.5

    # print(profHMM)
    transition_matrix = init_transition(profHMM, transition_matrix)
    pi = init_pi(all_states)
    # transition_matrix = convert_to_log(transition_matrix)
    # emission_matrix = convert_to_log(emission_matrix)
    # printing(transition_matrix)
    # printing(emission_matrix)
    # Training
    transition_matrix, emission_matrix = baum_welch(sequences_list, transition_matrix, emission_matrix, pi, is_DNA, 5)

    return transition_matrix, emission_matrix


def forward(observation, A, B, pi, residue_list):
    T = len(observation)
    c = np.zeros((T))
    Alpha = np.zeros((A.shape[0], T))
    Alpha[:, 0] = pi * B[:, residue_list.index(observation[0])]

    c[0] = 1.0 / np.sum(Alpha[:, 0])
    Alpha[:, 0] = c[0] * Alpha[:, 0]

    for t in range(1, T):  # time step
        Alpha[:, t] = np.dot(Alpha[:, t - 1], A) * B[:, residue_list.index(observation[t])]
        c[t] = 1.0 / np.sum(Alpha[:, t])
        Alpha[:, t] = Alpha[:, t] * c[t]

    log_Prob = -(np.sum(np.log(c)))

    return (log_Prob, Alpha, c)


def backward(observation, A, B, residue_list, c):
    T = len(observation)
    Beta = np.zeros((A.shape[0], T))
    Beta[:, T - 1] = 1.0
    Beta[:, T - 1] = Beta[:, T - 1] * c[T - 1]

    # backward from T-1 to 0
    for t in range(T - 2, -1, -1):
        Beta[:, t] = np.dot(A, (B[:, residue_list.index(observation[t + 1])] * Beta[:, t + 1]))
        Beta[:, t] = Beta[:, t] * c[t]
    return Beta


def baum_welch(observation_list, A, B, pi, is_DNA, n_iter):
    residue_list = DNA_list if is_DNA else amino_list
    states_num = A.shape[0]

    start_character = '<'
    end_character = '>'

    residue_list += [start_character, end_character]

    # to emit < and > at start and end
    newB = np.zeros((B.shape[0], B.shape[1] + 2))
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            newB[i, j] = B[i, j]
        if i == 0:
            newB[i, B.shape[1]] = 1
        elif i == B.shape[0] - 1:
            newB[i, B.shape[1] + 1] = 1
    B = newB

    M = B.shape[1]

    # to have pseudo-counts later/ we cannot say a probability is zero because it was not in the training set
    non_zeros_A = {}
    for state in range(A.shape[0]):
        non_zeros_A[state] = [s for s in range(A.shape[1]) if A[state, s] != 0]

    non_zeros_B = {}
    for state in range(B.shape[0]):
        non_zeros_B[state] = [s for s in range(B.shape[1]) if B[state, s] != 0]

    for n in range(n_iter):
        Expect_si_all = np.zeros(states_num)  # being in state i over all seqs
        Expect_si_all_TM1 = np.zeros(states_num)  # being in state i over all seqs until T-1
        Expect_si_sj_all_TM1 = np.zeros((states_num, states_num))
        Expect_si_t0_all = np.zeros(states_num)  # initially being in state i over all seqs
        Expect_si_vk_all = np.zeros((states_num, M))  # being in state i and seeing symbol vk

        for observation in observation_list:
            observation = list(observation)
            observation = [start_character] + observation + [end_character]
            log_Prob_Obs, Alpha, c = forward(observation, A, B, pi, residue_list)
            Beta = backward(observation, A, B, residue_list, c)
            T = len(observation)
            # Gamma = N * T matrix
            Gamma_raw = Alpha * Beta
            #  Gamma_raw.sum(0) = column summation (loop on states)
            Gamma = Gamma_raw / Gamma_raw.sum(0)
            Expect_si_t0_all += Gamma[:, 0]
            # Expect_si_all[ i ] = expected number of transitions from state i over all training sequences
            # Gamma.sum(1) = row summation (loop on Ts) for the denominator of emission matrix formula
            # For multiple observations: +=
            Expect_si_all += Gamma.sum(1)
            # Gamma[:, :T - 1].sum(1) = loop on Ts except the last for the denominator of transition matrix formula
            Expect_si_all_TM1 += Gamma[:, :T - 1].sum(1)
            # epsilon[ i,j,t ] = P(q_t = S_i, q_t+1 = S_j|Obs, hmm ) = N * N * T matrix
            epsilon = np.zeros((states_num, states_num, T - 1))

            for t in range(T - 1):
                for i in range(states_num):
                    epsilon[i, :, t] = Alpha[i, t] * A[i, :] * B[:, residue_list.index(observation[t+1])] * Beta[:, t+1]

            # for the nominator of transition matrix formula
            Expect_si_sj_all_TM1 += epsilon[:, :, :T - 1].sum(2)

            # for the nominator of emission matrix formula (only the observed characters)
            B_bar = np.zeros((states_num, M))
            for k in range(M):
                which = np.array([residue_list[k] == x for x in observation])
                B_bar[:, k] = Gamma.T[which, :].sum(0)
            # For multiple observations: +=
            Expect_si_vk_all += B_bar

        # Re-estimate
        Expect_si_t0_all = Expect_si_t0_all / np.sum(Expect_si_t0_all)
        pi = Expect_si_t0_all

        # transition
        A_bar = np.zeros((states_num ,states_num))
        for i in range(states_num):
            if Expect_si_all_TM1[i] == 0:
                A_bar[i, :] = A[i, :]
            else:
                A_bar[i, :] = Expect_si_sj_all_TM1[i, :] / Expect_si_all_TM1[i]
        A = A_bar

        # emission
        for i in range(states_num):
            if Expect_si_all[i] == 0:
                Expect_si_vk_all[i, :] = B[i, :]
            else:
                Expect_si_vk_all[i, :] = Expect_si_vk_all[i, :] / Expect_si_all[i]

        B = Expect_si_vk_all

        for state in range(A.shape[0]):
            for s in non_zeros_A[state]:
                if A[state, s] == 0:
                    A[state, s] = pseudo_count

        for state in range(B.shape[0]):
            for s in non_zeros_B[state]:
                if B[state, s] == 0:
                    B[state, s] = pseudo_count
    return A, B


def convert(s):
    new = ""
    for x in s:
        new += x
    return new


def viterbi(A, B, test_cases, is_DNA):
    state_count = A.shape[0]
    residue_list = DNA_list if is_DNA else amino_list

    gap_index_in_emission_matrix = residue_list.index('-')

    pruned_size = 4
    max_size = 8
    transitions_states = {}
    for state1 in range(state_count):
        transitions = [state2 for state2 in range(state_count) if A[state1, state2] != 0]
        transitions_states[state1] = transitions

    # start = time.time()
    for test_case in test_cases:
        test_case = "<" + test_case + ">"
        test_case = [char for char in test_case if char != '-']
        # path until here, number of characters up to here, prob of path up to here
        open_list = [[[0], 1, 1]]
        list_of_all_paths = []

        while len(open_list) > 0:
            if len(open_list) >= max_size:
                open_list.sort(key=lambda x: -x[2])
                open_list = open_list[0: pruned_size]

            current_node = open_list.pop()
            current_path = current_node[0]
            current_state = current_path[len(current_path) - 1]
            t = current_node[1]
            current_probability = current_node[2]

            if t == len(test_case):
                list_of_all_paths.append([current_path, current_probability])
                continue

            next_states = transitions_states[current_state]
            for state1 in next_states:
                # Delete
                if B[state1, gap_index_in_emission_matrix] == 1:
                    open_list.append([current_path + [state1], t, current_probability * A[current_state, state1]])
                else:
                    emission_probability = B[state1, residue_list.index(test_case[t])]
                    if emission_probability > 0:
                        open_list.append([current_path + [state1], t + 1,
                                          current_probability * A[current_state, state1] * emission_probability])

        list_of_all_paths.sort(key=lambda x: -x[1])

        # print("Number of paths found: ", len(closed_list))
        # print("Most likely path:")
        # print("************************************")
        # print()
        # for i in range(3):
        output_states = list_of_all_paths[0][0]
        test_case_counter = 1
        out_put = []
        for state_idx in range(1, len(output_states) - 1):
            if output_states[state_idx] % 3 == 0:  # delete
                out_put.append('-')
            elif output_states[state_idx] % 3 == 1:  # insert
                test_case_counter += 1
            else:  # match
                out_put.append(test_case[test_case_counter])
                test_case_counter += 1
        # print("Found Alignment: ", convert(out_put), "Probability: ", list_of_all_paths[i][1])
        # print("Elapsed time (Viterbi): ", time.time() - start, "seconds")
        return convert(out_put)


if __name__== "__main__":
    num_of_files, path = read_dir()
    for i in range(num_of_files):
        sequences_list, isDNA, test_cases, true_output = initialization(path, i)
        transition_matrix, emission_matrix = HMMStructure(sequences_list, isDNA)
        my_output = viterbi(transition_matrix, emission_matrix, test_cases, isDNA)
        insert_states = []
        sequences_bef_update = []
        DNA_list = DNA_list[0: 5]
        amino_list = amino_list[0: 21]
        if my_output == true_output:
            print(test_cases[0], ' ', true_output, ' ', my_output, ' True')
        else:
            print(test_cases[0], ' ', true_output, ' ', my_output, ' False')
