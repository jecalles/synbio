from copy import copy
from itertools import count
from math import prod
from queue import PriorityQueue
from collections import namedtuple

import numpy as np


# TODO: write tests
def get_overhangs(seq, L=4):
    # TODO: write docstring
    return set(
        seq[i:i + L] for i in range(len(seq) - L + 1)
    )


def find_overhangs(seq_list, L=4, min_solutions=1000):
    # TODO: write docstring

    item_count = count()
    Item = namedtuple("Item", "score item_num dist_matrix overhang_set".split())

    def _exit_condition(completed, ideal, min_solutions):
        best = completed.queue[0][2]
        return np.array_equal(best, ideal) or completed.qsize() >= min_solutions

    def _get_new_item(item, new_val, item_count):
        score, _, matrix, vals = copy(item)

        ix = len(vals)
        distances = [
            sum(char1 != char2 for char1, char2 in zip(new_val, val))
            for val in vals
        ]
        new_matrix = np.copy(matrix)
        new_matrix[:ix, ix] = distances

        new_score = -np.linalg.norm(new_matrix)
        new_vals = vals + [new_val]
        return Item(
            new_score, next(item_count), new_matrix, new_vals
        )

    # initialize some useful variables
    n = len(seq_list)
    overhangs_per_region = [get_overhangs(seq, L) for seq in seq_list]
    ideal = np.triu(np.ones((n, n)) * L, k=1)

    # handle min_solutions being larger than the solution space
    max_solutions = prod(
        len(overhang_set)
        for overhang_set in overhangs_per_region
    )
    num_solutions = min(min_solutions, max_solutions)

    # initialize queues for partial solutions and full solutions
    partial_solutions = PriorityQueue()
    completed_solutions = PriorityQueue()

    partial_solutions.put(
        Item(0, next(item_count), np.zeros((n, n)), [])
    )

    # Dijkstra's Algorithm
    while not partial_solutions.empty():
        item = partial_solutions.get()
        curr_region = len(item[-1])  # len of overhang set = ind of next region

        if curr_region == n:
            # base case: partial solution is complete
            completed_solutions.put(item)

            if _exit_condition(completed_solutions, ideal, num_solutions):
                return completed_solutions
        else:
            # recursive case: partial solution is incomplete
            new_items = (
                _get_new_item(item, new_val, item_count)
                for new_val in overhangs_per_region[curr_region]
                if new_val not in item[-1]
            )
            for it in new_items:
                partial_solutions.put(it)

    # we found no solutions
    raise RuntimeError("Unable to find appropriate overhang set(s)!")
