class Indel(object):

    def __init__(self, start_index, end_index, length, is_deletion):
        self.start_index = start_index
        self.end_index = end_index
        self.length = length
        self.is_deletion = is_deletion
        self.is_insertion = not is_deletion

    def envelope_cutsite(self, cutsite_index):
        return self.start_index <= cutsite_index and self.end_index >= cutsite_index
