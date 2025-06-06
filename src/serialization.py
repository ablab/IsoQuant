############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# serialization stuff
ENCODING = 'utf-8'
BYTE_ORDER = "big"
STR_LEN_BYTES = 4
NONE_STR_LEN = (1 << (8 * STR_LEN_BYTES)) - 1
SHORT_INT_BYTES = 2
LONG_INT_BYTES = 4
TERMINATION_INT = (1 << 32) - 1
SHORT_TERMINATION_INT = (1 << 16) - 1
SHORT_FLOAT_MULTIPLIER = 1 << 20
DICT_TYPE_LEN = 1
DICT_INT_TYPE = 9
DICT_INT_PAIR_TYPE = 10
DICT_STR_TYPE = 17


def write_string(s, outf):
    str_len = len(s)
    outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray(s, encoding=ENCODING))


def read_string(inf):
    str_len = int.from_bytes(inf.read(STR_LEN_BYTES), BYTE_ORDER)
    return inf.read(str_len).decode(encoding=ENCODING)


def write_string_or_none(s, outf):
    if s is None:
        outf.write(NONE_STR_LEN.to_bytes(STR_LEN_BYTES, BYTE_ORDER))
        return
    str_len = len(s)
    outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray(s, encoding=ENCODING))


def read_string_or_none(inf):
    str_len = int.from_bytes(inf.read(STR_LEN_BYTES), BYTE_ORDER)
    if str_len == NONE_STR_LEN:
        return None
    return inf.read(str_len).decode(encoding=ENCODING)


def write_int(val, outf, bytes_len=LONG_INT_BYTES):
    outf.write(val.to_bytes(bytes_len, BYTE_ORDER))


def read_int(inf, bytes_len=LONG_INT_BYTES):
    return int.from_bytes(inf.read(bytes_len), BYTE_ORDER)


def write_short_int(val, outf):
    write_int(val, outf, SHORT_INT_BYTES)


def read_short_int(inf):
    return read_int(inf, SHORT_INT_BYTES)


def write_list(l, outf, func):
    write_int(len(l), outf)
    for val in l:
        func(val, outf)


def read_list(inf, func):
    result = []
    list_size = read_int(inf)
    for i in range(list_size):
        result.append(func(inf))
    return result


def write_list_of_pairs(l, outf, func):
    write_int(len(l), outf)
    for val in l:
        func(val[0], outf)
        func(val[1], outf)


def read_list_of_pairs(inf, func):
    result = []
    list_size = read_int(inf)
    for i in range(list_size):
        result.append((func(inf), func(inf)))
    return result


def write_bool_array(bool_arr, outf):
    assert len(bool_arr) <= 8
    byte_val = 0
    for i, val in enumerate(bool_arr):
        if val:
            byte_val |= 1 << i
    outf.write(byte_val.to_bytes(1, BYTE_ORDER))


def read_bool_array(inf, arr_size=1):
    byte_val = int.from_bytes(inf.read(1), BYTE_ORDER)
    bool_arr = []
    for i in range(arr_size):
        mask = 1 << i
        bool_arr.append(byte_val & mask != 0)
    return bool_arr


# transforms a signed integer into an unsigned integer while using the 31st bit as a flag to remember the original sign
def write_int_neg(val, outf):
    if val < 0:
        val = abs(val)
        assert val & 1 << 31 == 0
        val |= 1 << 31
    else:
        assert val & 1 << 31 == 0
    outf.write(val.to_bytes(LONG_INT_BYTES, BYTE_ORDER))


def read_int_neg(inf):
    val = int.from_bytes(inf.read(LONG_INT_BYTES), BYTE_ORDER)
    if val & 1 << 31 != 0:
        return - (val & ((1 << 31) - 1))
    return val


def write_dict(d, outf):
    # keys - only strings, values - string, ints, int pairs
    write_int(len(d), outf)
    for k in d.keys():
        write_string(k, outf)
        v = d[k]
        if isinstance(v, int):
            outf.write(DICT_INT_TYPE.to_bytes(DICT_TYPE_LEN, BYTE_ORDER))
            write_int_neg(v, outf)
        elif isinstance(v, str):
            outf.write(DICT_STR_TYPE.to_bytes(DICT_TYPE_LEN, BYTE_ORDER))
            write_string(v, outf)
        elif isinstance(v, tuple):
            outf.write(DICT_INT_PAIR_TYPE.to_bytes(DICT_TYPE_LEN, BYTE_ORDER))
            write_int_neg(v[0], outf)
            write_int_neg(v[1], outf)
        else:
            raise ValueError("Unsupported value type in dictionary serialization")


def read_dict(inf):
    d = {}
    d_len = read_int(inf)
    for i in range(d_len):
        k = read_string(inf)
        vtype = int.from_bytes(inf.read(DICT_TYPE_LEN), BYTE_ORDER)
        if vtype == DICT_INT_TYPE:
            d[k] = read_int(inf)
        elif vtype == DICT_STR_TYPE:
            d[k] = read_string(inf)
        elif vtype == DICT_INT_PAIR_TYPE:
            d[k] = (read_int(inf), read_int(inf))
        else:
            raise ValueError("Serialized dictionary contains unsupported value")
    return d
