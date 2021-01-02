import itertools


class Location:
    # TODO: write docstring
    '''
    '''

    def __init__(self, start, end, strand='FWD'):
        # check that start and stop are valid
        if not start <= end:
            raise ValueError(f"start position ({start}) comes after end position ({end})")

        # check that strand input is valid
        if strand not in {'FWD', 'REV'}:
            raise ValueError(f"input strand ({strand}) not FWD or REV")

        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, {self.end}, {self.strand})"

    @staticmethod
    def contains(outer_loc, inner_loc):
        '''
        A static method that compares two locations, determining if the inner_loc is completely contained by the outer location.
        '''
        return (outer_loc.start <= inner_loc.start) and (inner_loc.end <= outer_loc.end)

    @staticmethod
    def overlaps(loc1, loc2):
        '''
        # TODO: write docstring for this method
        '''
        # TODO: implement lol
        return (loc1.start <= loc2.start and loc1.end > loc2.start) or (loc2.start <= loc1.start and loc2.end > loc1.start)

    @ staticmethod
    def find_overlaps(locations):
        '''
        A static method that, given a list of Location objects, constructs a graph representing all Location objects that overlap with each other.

        E.g., suppose our Locations are the following:


        a: ----------
        b:     ---
        c: ----
        d:           ----------
        e:      ----------

        >>> a = Location(0, 10)
        >>> b = Location(4, 7)
        >>> c = Location(0, 4)
        >>> d = Location(10, 20)
        >>> e = Location(5, 15)

        >>> Location.find_overlaps([a, b, c, d, e])
        [
            [1, 1, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 1, 0, 0],
            [0, 0, 0, 1, 1],
            [1, 1, 0, 1, 1],
        ]
        '''
        # calculate overlaps
        loc_pairs = itertools.product(locations, repeat=2)
        bool_array = [
            Location.overlaps(loc1, loc2) for
            (loc1, loc2) in loc_pairs
        ]
        # reshape array to square matrix
        n = len(locations)
        bool_matrix = [
            bool_array[i*n: (i+1)*n] for i in range(n)
        ]
        return bool_matrix

    def to_slice(self):
        '''
        '''
        # TODO: write docstring
        return slice(self.start, self.end, 1)


class Part:
    # TODO: write docstring
    '''
    '''

    def __init__(self, name=None, kind=None, seq="", location=None, subparts=[], metadata={}):
        pass

        # initialize location if not given
        if location is None:
            location = Location(0, len(seq))

        # check that subparts don't extend beyond this Part
        location_bool_array = [
            Location.contains(location, location) for subpart in subparts
        ]

        if not all(location_bool_array):
            overflowing_subparts = list(itertools.compress(subparts, location_bool_array))
            raise ValueError(
                f"The following subparts are located outside the range of this part {overflowing_subparts}")

        self.name = name
        self.kind = kind
        # TODO: Figure out seq retrieval that is strand-conscious and subsequence conscious (REV strand needs reverse indexing AND bp inversion)
        self.seq = seq
        self.location = location
        self.subparts = subparts
        self.metadata = metadata
