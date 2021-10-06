from typing import List, Tuple
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from code.fragments import start_stop_fragments
from code.codon_group import Group, codon_group_frequency
from code.phylip_matrix import draw_phylip_matrix


def main():
    codon_freq_tables: List[Tuple[str, dict[str, float]]] = []
    dicodon_freq_tables: List[Tuple[str, dict[str, float]]] = []
    file_paths = [
        str(filename)
        for filename in (Path(__file__).parent / "data").glob("*.fasta")
    ]

    for filename in file_paths:
        seq_record: SeqRecord = SeqIO.read(filename, "fasta")
        seq_fragments: List[Seq] = start_stop_fragments(seq_record)
        codon_freq_tables.append((
            seq_record.id,
            codon_group_frequency(seq_fragments, Group.CODON)
        ))
        dicodon_freq_tables.append((
            seq_record.id,
            codon_group_frequency(seq_fragments, Group.DICODON)
        ))

    draw_phylip_matrix(codon_freq_tables)
    draw_phylip_matrix(dicodon_freq_tables)


if __name__ == '__main__':
    main()
