import os

import pytest

from src import serialization

filename = "ser_test_file"


@pytest.fixture(scope="class")
def setup_class(request):
    request.cls.filehandler = open(filename, "xb")
    yield
    os.remove(filename)


@pytest.fixture(autouse=True)
def run_after_each_test(request):
    request.cls.filehandler = open(filename, "r+b")
    yield
    request.cls.filehandler.close()


@pytest.mark.usefixtures("setup_class")
class TestSerialization:

    @pytest.mark.parametrize(
        "func_wr, func_rd, value",
        [
            (serialization.write_int, serialization.read_int, 1),
            (serialization.write_int, serialization.read_int, 0),
            (serialization.write_int, serialization.read_int, 65536),
            (serialization.write_short_int, serialization.read_short_int, 0),
            (serialization.write_short_int, serialization.read_short_int, 65535),
            (serialization.write_string, serialization.read_string, ""),
            (serialization.write_string, serialization.read_string, "quite a long string dont you think"),
            (serialization.write_string_or_none, serialization.read_string_or_none, None),
            (serialization.write_string_or_none, serialization.read_string_or_none, "None"),
            (serialization.write_int_neg, serialization.read_int_neg, -1),
            (serialization.write_int_neg, serialization.read_int_neg, -65536),
            (serialization.write_bool_array, serialization.read_bool_array, [False]),
            (serialization.write_dict, serialization.read_dict, {'key1': 'value1', 'key2': 'value2'}),
            (serialization.write_dict, serialization.read_dict, {'key1': 0, 'key2': (655536, 1)}),
            (serialization.write_dict, serialization.read_dict, {})
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
            (serialization.write_int, serialization.read_int, [65539, 0, 1]),
            (serialization.write_int, serialization.read_int, []),
            (serialization.write_string_or_none, serialization.read_string_or_none, ["str", "", None]),
            (serialization.write_string_or_none, serialization.read_string_or_none, []),
            (serialization.write_string, serialization.read_string, ["str", "quite a long string dont you think", ""]),
            (serialization.write_string, serialization.read_string, []),
            (serialization.write_int_neg, serialization.read_int_neg, [-6, -33333, 0]),
            (serialization.write_int_neg, serialization.read_int_neg, [])
        ]
    )
    def test_write_read_lists_positive(self, func_wr, func_rd, value):
        serialization.write_list(value, self.filehandler, func_wr)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = serialization.read_list(self.filehandler, func_rd)
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
        serialization.write_bool_array(value, self.filehandler)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = serialization.read_bool_array(self.filehandler, len(value))
        assert actual == value

    @pytest.mark.parametrize(
        "func_wr, func_rd, list_of_pairs",
        [
            (serialization.write_int, serialization.read_int, []),
            (serialization.write_int, serialization.read_int, [(1225, 78854)]),
            (serialization.write_int, serialization.read_int, [(1225, 78854), (1, 0), (65536, 65536)]),
            (serialization.write_string, serialization.read_string, [("1225", "78854"), ("", " "), ("65536", "65536")])
        ]
    )
    def test_write_read_list_of_pairs_positive(self, func_wr, func_rd, list_of_pairs):
        serialization.write_list_of_pairs(list_of_pairs, self.filehandler, func_wr)
        self.filehandler.flush()
        self.filehandler.seek(0)
        actual = serialization.read_list_of_pairs(self.filehandler, func_rd)
        assert actual == list_of_pairs
