from typing import List, Set
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from code.codon_group import Group, codon_group


def start_stop_fragments(
    seq_record: SeqRecord,
    minimum_fragment_length: int = 100,
    start_codons: Set[str] = {"ATG"},
    stop_codons: Set[str] = {"TAA", "TAG", "TGA"}
) -> List[Seq]:
    nucleotides: Seq = seq_record.seq
    start_stop_fragments: List[Seq] = []
    nucleotide_strand: Seq
    for nucleotide_strand in [nucleotides, nucleotides.reverse_complement()]:
        for frame in range(3):
            frame_codons: List[str] = codon_group(
                nucleotide_strand, Group.CODON, frame)
            fragment = Seq("")
            start_reading = False
            for codon in frame_codons:
                if codon in start_codons:
                    start_reading = True
                    fragment += codon
                elif start_reading:
                    fragment += codon
                    if codon in stop_codons:
                        start_reading = False
                        start_stop_fragments.append(fragment)
                        fragment = Seq("")

    return [
        fragment for fragment in start_stop_fragments
        if len(fragment) > minimum_fragment_length
    ]
