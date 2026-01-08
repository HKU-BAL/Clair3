import os
import collections


_Region = collections.namedtuple('Region', 'ref_name start end')

class Region(_Region):
    """Represents a genomic region."""

    @property
    def name(self):
        """Samtools-style region string, zero-base end exclusive."""
        return self.__str__()

    def __str__(self):
        """Return string representation of region."""
        # This will be zero-based, end exclusive
        start = 0 if self.start is None else self.start
        end = '' if self.end is None else self.end
        return '{}:{}-{}'.format(self.ref_name, start, end)

    @property
    def size(self):
        """Return size of region."""
        return self.end - self.start

    @classmethod
    def from_string(cls, region):
        """Parse region string into `Region` objects.

        :param region: region str

        >>> Region.from_string('Ecoli') == Region(
        ...     ref_name='Ecoli', start=None, end=None)
        True
        >>> Region.from_string('Ecoli:1000-2000') == Region(
        ...     ref_name='Ecoli', start=1000, end=2000)
        True
        >>> Region.from_string('Ecoli:1000') == Region(
        ...     ref_name='Ecoli', start=1000, end=None)
        True
        >>> Region.from_string('Ecoli:-1000') == Region(
        ...     ref_name='Ecoli', start=0, end=1000)
        True
        >>> Region.from_string('Ecoli:500-') == Region(
        ...     ref_name='Ecoli', start=500, end=None)
        True
        >>> Region.from_string('A:B:c:500-') == Region(
        ...     ref_name='A:B:c', start=500, end=None)
        True
        """
        if ':' not in region:
            ref_name, start, end = region, None, None
        else:
            start, end = None, None
            ref_name, bounds = region.rsplit(':', 1)
            if bounds[0] == '-':
                start = 0
                end = int(bounds.replace('-', ''))
            elif '-' not in bounds:
                start = int(bounds)
                end = None
            elif bounds[-1] == '-':
                start = int(bounds[:-1])
                end = None
            else:
                start, end = [int(b) for b in bounds.split('-')]
        return cls(ref_name, start, end)

    def split(region, size, overlap=0, fixed_size=True):
        """Split region into sub-regions of a given length.

        :param size: size of sub-regions.
        :param overlap: overlap between ends of sub-regions.
        :param fixed_size: ensure all sub-regions are equal in size. If `False`
            then the final chunk will be created as the smallest size to
            conform with `overlap`.

        :returns: a list of sub-regions.

        """
        regions = list()
        if size >= region.size:
            return [region]
        for start in range(region.start, region.end, size - overlap):
            end = min(start + size, region.end)
            regions.append(Region(region.ref_name, start, end))
        if len(regions) > 1:
            if fixed_size and regions[-1].size < size:
                del regions[-1]
                end = region.end
                start = end - size
                if start > regions[-1].start:
                    regions.append(Region(region.ref_name, start, end))
        return regions