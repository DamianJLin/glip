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

    def reverse(self, style='fix_0'):
        """
        Reverse elements in place. Style 'fix_0' fixes the element at position 0, and 'reflect'
        reflects.
        """

        if style == 'fix_0':
            self.elements.reverse()
            self.elements = self.elements[-1:] + self.elements[:-1]
        elif style == 'reflect':
            self.elements.reverse()
        else:
            raise ValueError("Style must be 'fix_0' or 'reflect'")

    def cycle(self, offset=1):
        self.elements = self.elements[offset:] + self.elements[:offset]
