############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
import io
from src import serialization
from src.serialization import *


@pytest.fixture
def setup_class(request):
    request.cls.filehandler = io.BytesIO()

    def teardown():
        request.cls.filehandler.close()

    request.addfinalizer(teardown)


@pytest.mark.usefixtures("setup_class")
class TestSerialization:
    @pytest.mark.parametrize(
        "func_wr, func_rd, value",
        [
            (write_int, read_int, 1),
            (write_int, read_int, 0),
            (write_int, read_int, 65536),
            (write_short_int, read_short_int, 0),
            (write_short_int, read_short_int, 65535),
            (write_string, read_string, ""),
            (write_string, read_string, "quite a long string dont you think"),
            (write_string_or_none, read_string_or_none, None),
            (write_string_or_none, read_string_or_none, "None"),
            (write_int_neg, read_int_neg, -1),
            (write_int_neg, read_int_neg, -65536),
            (write_bool_array, read_bool_array, [False]),
            (write_dict, read_dict, {'key1': 'value1', 'key2': 'value2'}),
            (write_dict, read_dict, {'key1': 0, 'key2': (655536, 1)}),
            (write_dict, read_dict, {})
        ]
    )
    def test_write_read_positive(self, value, func_wr, func_rd):
        func_wr(value, self.filehandler)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = func_rd(self.filehandler)
        assert actual == value

    @pytest.mark.parametrize(
        "func_wr, func_rd, value",
        [
            (write_int, read_int, [65539, 0, 1]),
            (write_int, read_int, []),
            (write_string_or_none, read_string_or_none, ["str", "", None]),
            (write_string_or_none, read_string_or_none, []),
            (write_string, read_string, ["str", "quite a long string dont you think", ""]),
            (write_string, read_string, []),
            (write_int_neg, read_int_neg, [-6, -33333, 0]),
            (write_int_neg, read_int_neg, [])
        ]
    )
    def test_write_read_lists_positive(self, func_wr, func_rd, value):
        write_list(value, self.filehandler, func_wr)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = read_list(self.filehandler, func_rd)
        assert actual == value

    @pytest.mark.parametrize(
        "value",
        [
            ([]),
            ([True, False]),
            ([False, False, False, True, False, True, True, False]),
        ]
    )
    def test_write_read_bool_array_positive(self, value):
        write_bool_array(value, self.filehandler)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = read_bool_array(self.filehandler, len(value))
        assert actual == value

    @pytest.mark.parametrize(
        "func_wr, func_rd, list_of_pairs",
        [
            (write_int, read_int, []),
            (write_int, read_int, [(1225, 78854)]),
            (write_int, read_int, [(1225, 78854), (1, 0), (65536, 65536)]),
            (write_string, read_string, [("1225", "78854"), ("", " "), ("65536", "65536")])
        ]
    )
    def test_write_read_list_of_pairs_positive(self, func_wr, func_rd, list_of_pairs):
        write_list_of_pairs(list_of_pairs, self.filehandler, func_wr)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = read_list_of_pairs(self.filehandler, func_rd)
        assert actual == list_of_pairs

    def test_read_write_list_with_custom_serializer(self):
        class TestClass:
            def __init__(self, value):
                self.value = value

            def serialize(self, outfile):
                write_int(self.value, outfile)

            @staticmethod
            def deserialize(infile):
                obj = TestClass(0)
                obj.value = read_int(infile)
                return obj

        test_objects = [TestClass(1), TestClass(2), TestClass(3)]
        write_list(test_objects, self.filehandler, TestClass.serialize)
        self.filehandler.flush()
        self.filehandler.seek(0)
        read_objects = read_list(self.filehandler, TestClass.deserialize)

        assert len(test_objects) == len(read_objects)
        for orig, read in zip(test_objects, read_objects):
            assert orig.value == read.value

    def test_termination_int(self):
        write_int(TERMINATION_INT, self.filehandler)
        self.filehandler.flush()
        self.filehandler.seek(0)
        value = read_int(self.filehandler)
        assert value == TERMINATION_INT

