from indel_mapper_lib.indel import Indel

class Read(object):

    def __init__(self, query_name, reference_name, reference_positions, query_sequence, aligned_pairs):
        self.query_name = query_name
        self.query_sequence = query_sequence
        self.aligned_pairs = aligned_pairs
        self.reference_name = reference_name
        self.reference_positions = reference_positions
        self.indels = self._compute_indels()
        self.valid_indels = self._valid_indels()

    def _compute_indels(self):
        indels = []

        start_indel_index = self.reference_positions[0]
        prev_reference_index = self.reference_positions[0]
        indel_length = 0

        for idx in self.reference_positions[1:]:
            if idx is not None:
                if prev_reference_index is not None:
                    # deletion
                    if idx != prev_reference_index + 1 and idx != prev_reference_index:
                        indels.append(Indel(prev_reference_index, idx, (idx-1) - prev_reference_index, True))
                else:
                    # end of insertion
                    indels.append(Indel(start_indel_index, idx, indel_length, False))
                    indel_length = 0
            else:
                # start of insertion
                if prev_reference_index is not None:
                    start_indel_index = prev_reference_index
                    indel_length += 1
                else:
                    indel_length += 1

            prev_reference_index = idx

        return indels

    def _valid_indels(self):
        return [indel for indel in self.indels if (indel.start_index is not None and indel.end_index is not None)]
