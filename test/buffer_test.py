# -*- coding: UTF-8 -*-
import unittest
from lib.buffer import *


class TestArrayBuffer(unittest.TestCase):
    def test_setitem_getitem(self):
        p = ArrayBuffer(5, 5, -123)
        p[0] = 0
        p[1] = 10
        p[3] = 30
        p[2] = 20
        p[4] = 40
        p[0] += 1
        p[1] += 2
        p[2] += 3
        p[3] += 4
        p[4] += 5
        self.assertEqual(p[0], 1)
        self.assertEqual(p[1], 12)
        self.assertEqual(p[2], 23)
        self.assertEqual(p[3], 34)
        self.assertEqual(p[4], 45)

    def test_pop_one(self):
        p = ArrayBuffer(10, 5, -123)
        values = [x for x in p.pop_one()]
        self.assertEqual(values, [-123])
        self.assertEqual(p.pos, 1)

        values = [x for x in p.pop_one()]
        self.assertEqual(values, [-123])
        self.assertEqual(p.pos, 2)

        p[2] = 20
        p[5] = 50
        p[3] = 30
        p[4] = 40
        p[2] = 21
        p[3] = 31

        values = [x for x in p.pop_one()]
        self.assertEqual(values, [21])
        self.assertEqual(p.pos, 3)

        values = [x for x in p.pop_one()]
        self.assertEqual(values, [31])
        self.assertEqual(p.pos, 4)

        values = [x for x in p.pop_one()]
        self.assertEqual(values, [40])
        self.assertEqual(p.pos, 5)

        self.assertEqual([x for x in p.pop_one()], [50])
        self.assertEqual([x for x in p.pop_one()], [-123])
        self.assertEqual([x for x in p.pop_one()], [-123])
        self.assertEqual([x for x in p.pop_one()], [-123])

        p[9] = 90
        values = [x for x in p.pop_one()]
        self.assertEqual(values, [90])
        self.assertEqual(p.pos, 10)

        values = [x for x in p.pop_one()]
        self.assertEqual(values, [])
        self.assertEqual(p.pos, 10)

    def test_pop_until(self):
        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        values = [x for x in p.pop_until(5)]
        values_bis = [x for x in p.pop_until(5)]
        self.assertEqual(values, [0, 1, 2, 3, 4])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 5)

        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_until(10)]
        values_bis = [x for x in p.pop_until(10)]
        self.assertEqual(values, [5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_until(5)]
        values_bis = [x for x in p.pop_until(5)]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 5)

        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        for _ in range(10): list(p.pop_one())
        values = [x for x in p.pop_until(10)]
        values_bis = [x for x in p.pop_until(10)]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = ArrayBuffer(0, 10)
        values = [x for x in p.pop_until(0)]
        values_bis = [x for x in p.pop_until(0)]
        self.assertEqual(values, [])
        self.assertEqual(p.pos, 0)

    def test_pop_all(self):
        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = ArrayBuffer(10, 10)
        for i in range(10): p[i] = i
        for _ in range(10): list(p.pop_one())
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = ArrayBuffer(0, 10)
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [])
        self.assertEqual(p.pos, 0)

        p = ArrayBuffer(10, 5, -123)
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [-123] * 10)
        self.assertEqual(values_bis, [])


class TestPrefixSumBuffer(unittest.TestCase):
    def test_pop_one(self):
        p = PrefixSumBuffer(10, 5, 123)
        values = [x for x in p.pop_one()]
        values_bis = [x for x in p.pop_one()]
        self.assertEqual(values, [123])
        self.assertEqual(values_bis, [123])

        p = PrefixSumBuffer(10, 5)
        values = [x for x in p.pop_one()]
        self.assertEqual(values, [0])

        p = PrefixSumBuffer(1, 1)
        values = [x for x in p.pop_one()]
        values_bis = [x for x in p.pop_one()]
        self.assertEqual(values, [0])
        self.assertEqual(values_bis, [])

    def test_add_interval(self):
        p = PrefixSumBuffer(10, 10)
        p.add_interval(0, 10, 5)
        values = [x for x in p.pop_all()]
        self.assertEqual(values, [5] * 10)

        p = PrefixSumBuffer(10, 10, 3)
        p.add_interval(5, 7, 10)
        values = [x for x in p.pop_all()]
        self.assertEqual(values, [3, 3, 3, 3, 3, 13, 13, 3, 3, 3])

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        values = [x for x in p.pop_all()]
        self.assertEqual(values, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    def test_pop_until(self):
        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        values = [x for x in p.pop_until(5)]
        values_bis = [x for x in p.pop_until(5)]
        self.assertEqual(values, [0, 1, 2, 3, 4])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 5)

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_until(10)]
        values_bis = [x for x in p.pop_until(10)]
        self.assertEqual(values, [5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_until(5)]
        values_bis = [x for x in p.pop_until(5)]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 5)

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        for _ in range(10): list(p.pop_one())
        values = [x for x in p.pop_until(10)]
        values_bis = [x for x in p.pop_until(10)]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = PrefixSumBuffer(0, 10)
        values = [x for x in p.pop_until(0)]
        values_bis = [x for x in p.pop_until(0)]
        self.assertEqual(values, [])
        self.assertEqual(p.pos, 0)

    def test_pop_all(self):
        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        for _ in range(5): list(p.pop_one())
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [5, 6, 7, 8, 9])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = PrefixSumBuffer(10, 10)
        for i in range(10): p.add_interval(i+1, 10)
        for _ in range(10): list(p.pop_one())
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [])
        self.assertEqual(values_bis, [])
        self.assertEqual(p.pos, 10)

        p = PrefixSumBuffer(0, 10)
        values = [x for x in p.pop_all()]
        values_bis = [x for x in p.pop_all()]
        self.assertEqual(values, [])
        self.assertEqual(p.pos, 0)

        p = PrefixSumBuffer(10, 5)
        values = [x for x in p.pop_all()]
        self.assertEqual(values, [0] * 10)
