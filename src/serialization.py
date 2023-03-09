# serialization stuff
ENCODING = 'utf-8'
BYTE_ORDER = "big"
STR_LEN_BYTES = 2
SHORT_INT_BYTES = 2
LONG_INT_BYTES = 3


def write_string(s, outf):
    str_len = len(s)
    outf.write(str_len.to_bytes(STR_LEN_BYTES, BYTE_ORDER) + bytearray(s, encoding=ENCODING))


def read_string(inf):
    str_len = int.from_bytes(inf.read(STR_LEN_BYTES), BYTE_ORDER)
    return inf.read(str_len).decode(encoding=ENCODING)


def write_int(val, outf, bytes_len=LONG_INT_BYTES):
    outf.write(val.to_bytes(bytes_len, BYTE_ORDER))


def read_int(inf, bytes_len=LONG_INT_BYTES):
    return int.from_bytes(inf.read(bytes_len), BYTE_ORDER)


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



