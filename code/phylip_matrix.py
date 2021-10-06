from typing import List, Tuple
from math import sqrt
from statistics import mean


def calculate_distance_between(
    one: dict[str, float],
    another: dict[str, float]
) -> float:
    return sqrt(mean([
        pow(one[key] - another[key], 2)
        for key in one
    ]))


def draw_phylip_matrix(
    freq_tables: List[Tuple[str, dict[str, float]]]
):
    print(f'\n{len(freq_tables)}')

    for freq_table in freq_tables:
        print(freq_table[0], end=" ")
        for other_freq_table in freq_tables:
            distance = calculate_distance_between(
                freq_table[1],
                other_freq_table[1]
            )
            print(f'{format(distance, ".6f")}', end=" ")
        print("")
    print("")
