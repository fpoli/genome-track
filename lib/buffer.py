# -*- coding: UTF-8 -*-


class ArrayBuffer:
    """
    Array-like data structure that allows set/get operations only on elements
    whose index is in a sliding buffer that starts at [0, buffer_size - 1] and
    slowly moves to the end of the array.

    Ranges:
        array range: [0, size - 1]
        initial buffer range: [0, buffer_size - 1]
        generic buffer range: [pos, pos + buffer_size - 1]
        final buffer range:   [size - buffer_size, size - 1]
        ... then buffer_size is decreased until it reaches 0
    """
    def __init__(self, size, buffer_size, default=0):
        self.size = size
        self.buffer_size = buffer_size
        self.default = default
        self.buffer_start = 0
        self.buffer_end = self.buffer_size
        self.circular_buffer = [self.default] * self.buffer_size

    @property
    def pos(self):
        return self.buffer_start

    def __getitem__(self, pos):
        """Return the value of the element in position pos, if it is in the
        buffer range.
        """
        assert self.buffer_start <= pos < self.buffer_end
        return self.circular_buffer[pos % self.buffer_size]

    def __setitem__(self, pos, value):
        """Set the value of the element in position pos, if it is in the
        buffer range.
        """
        assert self.buffer_start <= pos < self.buffer_end
        self.circular_buffer[pos % self.buffer_size] = value

    def pop_one(self):
        """Generate the element with lowest index in the buffer, then increase
        the buffer range edges by one position.
        """
        assert self.buffer_start <= self.buffer_end <= self.size

        if self.buffer_start < self.buffer_end:
            yield self[self.buffer_start]
            self.buffer_start += 1

            if self.buffer_end < self.size:
                self.buffer_end += 1
                self[self.buffer_end - 1] = self.default

        assert self.buffer_start <= self.buffer_end <= self.size

    def pop_until(self, pos):
        """Generate elements while moving the buffer to index 'pos'.
        """
        assert self.buffer_start <= pos <= self.size

        while self.buffer_start < pos:
            yield self[self.buffer_start]
            self.buffer_start += 1

            if self.buffer_end < self.size:
                self.buffer_end += 1
                self[self.buffer_end - 1] = self.default

        assert self.buffer_start == pos

    def pop_all(self):
        """Generate elements while moving the buffer to the end of the array.
        """
        return self.pop_until(self.size)


class PrefixSumBuffer:
    """
    Array-like data structure that allows add_interval operations only on
    elements whose index is in a sliding buffer that starts at
    [0, buffer_size - 1] and slowly moves to the end of the array.

    Ranges:
        array range: [0, size - 1]
        initial buffer range: [0, buffer_size - 1]
        generic buffer range: [pos, pos + buffer_size - 1]
        last buffer range:   [size - buffer_size, size - 1]
        ... then buffer_size is decreased until it reaches 0
    """
    def __init__(self, size, buffer_size, default=0):
        self.size = size
        self.buffer = ArrayBuffer(size, min(size, buffer_size + 1))
        if size > 0:
            self.buffer[0] = default
        self.current = 0

    @property
    def pos(self):
        return self.buffer.pos

    def add_interval(self, start, end, value=1):
        """Increment elements with index between start and end - 1 (included).
        """
        assert start <= end
        assert self.buffer.buffer_start <= start
        assert end <= self.buffer.buffer_start + self.buffer.buffer_size

        if start < self.buffer.size:
            self.buffer[start] += value
            if end < self.buffer.size:
                self.buffer[end] -= value

    def pop_one(self):
        """Generate the element with lowest index in the buffer, then increase
        the buffer range edges by one position.
        """
        for value in self.buffer.pop_one():
            self.current += value
            yield self.current

    def pop_until(self, pos):
        """Generate elements while moving the buffer to index 'pos'.
        """
        assert self.buffer.buffer_start <= pos <= self.size

        for value in self.buffer.pop_until(pos):
            self.current += value
            yield self.current

        assert self.buffer.buffer_start == pos

    def pop_all(self):
        """Generate elements while moving the buffer to the end of the array.
        """
        return self.pop_until(self.size)
