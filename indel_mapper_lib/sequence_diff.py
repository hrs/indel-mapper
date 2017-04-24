from enum import Enum

from .diff_components import DiffDeletion
from .diff_components import DiffInsertion
from .diff_components import DiffMatch
from .diff_components import DiffMetadata


class DiffState(Enum):
    MATCH = 1
    INSERTION = 2
    DELETION = 3
    METADATA = 4


class SequenceDiff(object):
    STATE_MAP = {
        DiffState.DELETION: DiffDeletion,
        DiffState.INSERTION: DiffInsertion,
        DiffState.MATCH: DiffMatch,
        DiffState.METADATA: DiffMetadata,
    }

    def __init__(self, original, changed):
        self._validate(original, changed)
        self.changes = self._change_tokens(original, changed)

    def description(self):
        relevant_changes = self._strip_end_matches(self.changes)
        desc = ""

        for first, second in zip(relevant_changes, relevant_changes[1:]):
            desc += str(first)
            if first.is_metadata() or second.is_metadata():
                desc += " "
            else:
                desc += ", "

        desc += str(relevant_changes[-1])

        return desc

    def _change_tokens(self, original, changed):
        curr = ""
        state = DiffState.MATCH
        tokens = []

        # Compare each pair of bases
        for base_1, base_2 in zip(original, changed):
            # Begin/continue reading a DiffMetadata
            if base_1 == "|" and base_2 == "|":
                if state == DiffState.METADATA:
                    curr += "|"
                else:
                    tokens.append(self._create_diff(curr, state))
                    curr = "|"
                    state = DiffState.METADATA

            # Begin/continue reading a DiffMatch
            elif base_1 == base_2:
                if state == DiffState.MATCH:
                    curr += base_1
                else:
                    tokens.append(self._create_diff(curr, state))
                    curr = base_1
                    state = DiffState.MATCH

            # Begin/continue reading a DiffInsertion
            elif base_1 == "_":
                if state == DiffState.INSERTION:
                    curr += base_2
                else:
                    tokens.append(self._create_diff(curr, state))
                    curr = base_2
                    state = DiffState.INSERTION

            # Begin/continue reading a DiffDeletion
            elif base_2 == "_":
                if state == DiffState.DELETION:
                    curr += base_1
                else:
                    tokens.append(self._create_diff(curr, state))
                    curr = base_1
                    state = DiffState.DELETION

            else:
                raise ValueError("Can't diff bases \"{}\" and \"{}\"".format(base_1, base_2))

        # Append the last token
        if curr != "":
            tokens.append(self._create_diff(curr, state))

        return tokens

    def _create_diff(self, sequence, state):
        if state in self.STATE_MAP:
            return self.STATE_MAP[state](sequence)
        else:
            raise ValueError("Unknown diff state \"{}\"".format(state))

    def _validate(self, original, changed):
        if len(original) != len(changed):
            raise ValueError("Can't diff sequences of differing length")
        if self._cutsite_positions(original) != self._cutsite_positions(changed):
            raise ValueError("Can't diff sequences with misaligned cutsites")

    def _cutsite_positions(self, seq):
        return [pos for pos, base in enumerate(seq) if base == "|"]

    def _strip_end_matches(self, changes):
        while len(changes) != 0 and changes[0].is_match():
            changes = changes[1:]
        while len(changes) != 0 and changes[-1].is_match():
            changes = changes[:-1]
        return changes
