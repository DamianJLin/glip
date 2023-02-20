class CyclicList():

    def __init__(self, *elements):

        self.elements = list(elements)
        self.n = len(self.elements)

    def __getitem__(self, key):

        if not isinstance(key, int):
            raise TypeError(
                f'CyclicList indexes must be int, not {type(key)}.'
            )

        return self.elements[key % self.n]

    def __setitem__(self, key, value):

        if not isinstance(key, int):
            raise TypeError(
                f'CyclicList indexes must be int, not {type(key)}.'
            )

        self.elements[key % self.n] = value

    def __repr__(self):

        return '-' + str(self.elements) + '-'

    def append(self, element):

        self.elements.append(element)
        self.n = len(self.elements)

    def __len__(self):

        return self.n

    def __iter__(self):

        yield from self.elements
