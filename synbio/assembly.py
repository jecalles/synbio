from typing import List, Iterator
from dataclasses import dataclass
from itertools import count
from math import prod
from queue import PriorityQueue

import numpy as np

from synbio.interfaces import SeqType
from synbio.utils import reverse_complement


@dataclass(order=True)
class OverhangSet:
    score: int
    item_num: int
    dist_matrix: np.ndarray
    overhangs: List[SeqType]


# TODO: write tests
def get_overhangs(seq, L=4):
    # TODO: write docstring
    all_overhangs = set(
        str(seq[i:i + L])
        for i in range(len(seq) - L + 1)
    )
    minus_self_complimentary = set(
        overhang for overhang in all_overhangs
        if overhang != reverse_complement(overhang)
    )
    return minus_self_complimentary


def find_overhangs(
        seq_list: List[SeqType],
        L: int = 4,
        min_solutions: int = 1000
) -> List[OverhangSet]:
    def _exit_condition(
            completed: PriorityQueue,
            ideal: int,
            min_solutions: int
    ) -> bool:
        best_solution = completed.queue[0]
        if not isinstance(best_solution, OverhangSet):
            raise TypeError
        best = best_solution.dist_matrix

        return np.array_equal(best, ideal) or completed.qsize() >= min_solutions

    def _get_new_item(
            item: OverhangSet,
            new_val: SeqType,
            counter: Iterator[int]
    ) -> OverhangSet:
        matrix_ix = len(item.overhangs)
        distances = [
            sum(char1 != char2 for char1, char2 in zip(new_val, val))
            for val in item.overhangs
        ]
        new_matrix = np.copy(item.dist_matrix)
        new_matrix[:matrix_ix, matrix_ix] = distances

        new_score = -np.linalg.norm(new_matrix) # TODO: update scoring alg.
        new_vals = item.overhangs + [new_val]
        return OverhangSet(
            new_score, next(counter), new_matrix, new_vals
        )

    # initialize some useful variables
    item_count = count()
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
        OverhangSet(0, next(item_count), np.zeros((n, n)), [])
    )

    # Dijkstra's Algorithm
    while not partial_solutions.empty():
        item = partial_solutions.get()
        curr_region = len(item.overhangs)  # len of overhang set = ind of
        # next region

        if curr_region == n:
            # base case: partial solution is complete
            completed_solutions.put(item)

            if _exit_condition(completed_solutions, ideal, num_solutions):
                return completed_solutions.queue
        else:
            # recursive case: partial solution is incomplete
            new_items = (
                _get_new_item(item, new_val, item_count)
                for ix, new_val in enumerate(overhangs_per_region[curr_region])
                if new_val not in item.overhangs
            )
            for it in new_items:
                partial_solutions.put(it)

    # we found no solutions
    raise RuntimeError("Unable to find appropriate overhang set(s)!")
