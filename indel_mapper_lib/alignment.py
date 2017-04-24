class Alignment(object):
    def __init__(self, reference, read):
        self.reference = reference
        self.read = read

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (other.reference == self.reference and
                    other.read == self.read)
        return False

    def __repr__(self):
        return "({}, {})".format(self.reference, self.read)
