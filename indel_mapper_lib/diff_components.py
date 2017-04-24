class DiffComponent(object):
    def __init__(self, sequence):
        self.sequence = sequence

    def is_match(self):
        return False

    def is_metadata(self):
        return False

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                other.sequence == self.sequence)


class DiffMatch(DiffComponent):
    def is_match(self):
        return True

    def __repr__(self):
        name = "match"
        if len(self.sequence) > 1:
            name = "matches"
        return "{} {}".format(len(self.sequence), name)


class DiffInsertion(DiffComponent):
    def __repr__(self):
        name = "insertion"
        if len(self.sequence) > 1:
            name = "insertions"
        return "{} {} ({})".format(len(self.sequence), name, self.sequence)


class DiffDeletion(DiffComponent):
    def __repr__(self):
        name = "deletion"
        if len(self.sequence) > 1:
            name = "deletions"
        return "{} {} ({})".format(len(self.sequence), name, self.sequence)


class DiffMetadata(DiffComponent):
    def is_metadata(self):
        return True

    def __repr__(self):
        return self.sequence
