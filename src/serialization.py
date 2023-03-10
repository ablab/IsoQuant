# serialization stuff
ENCODING = 'utf-8'
BYTE_ORDER = "big"
STR_LEN_BYTES = 2
SHORT_INT_BYTES = 2
LONG_INT_BYTES = 4
TERMINATION_INT = (1 << 32) - 1
SHORT_TERMINATION_INT = (1 << 16) - 1
SHORT_FLOAT_MULTIPLIER = 1 << 20


def write_string(s, outf):
    str_len = len(s)
    outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray(s, encoding=ENCODING))


def read_string(inf):
    str_len = int.from_bytes(inf.read(STR_LEN_BYTES), BYTE_ORDER)
    return inf.read(str_len).decode(encoding=ENCODING)


def write_string_or_none(s, outf):
    if s is None:
        str_len = 1
        outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray([0]))
        return
    str_len = len(s)
    outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray(s, encoding=ENCODING))


def read_string_or_none(inf):
    str_len = int.from_bytes(inf.read(STR_LEN_BYTES), BYTE_ORDER)
    if str_len == 1:
        val = inf.read(str_len)
        if val[0] == 0:
            return None
        else:
            return val.decode(encoding=ENCODING)
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