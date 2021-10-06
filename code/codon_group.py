from enum import Enum
from itertools import combinations_with_replacement, permutations
from typing import Dict, List
from Bio.Seq import Seq


class Group(Enum):
    CODON = 1
    DICODON = 2


def codon_group(
    seq: Seq,
    group: Group = Group.CODON,
    frame: int = 0
) -> List[str]:
    length: int = 3 * group.value
    return [
        str(chunk)
        for index in range(frame, len(seq), length)
        if len(chunk := seq[index:index + length]) == length
    ]


def codon_group_frequency(
    fragments: List[Seq],
    group: Group,
    nucleobases: str = "ATGC"
) -> dict[str, float]:
    possible_codons: List[str] = list({
        ''.join(permutation)
        for combination in combinations_with_replacement(
            nucleobases, 3 * group.value)
        for permutation in permutations(combination)
    })
    frequency_dict: Dict[str, int] = {
        codon: 0
        for codon in possible_codons
    }

    total_codons: int = 0
    for frag in fragments:
        frag_codons = codon_group(frag, group)
        total_codons += len(frag_codons)
        for codon in possible_codons:
            frequency_dict[codon] += frag_codons.count(codon)

    return {
        key: (value / total_codons)
        for key, value in frequency_dict.items()
    }
