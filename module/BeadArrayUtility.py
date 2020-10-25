import struct
from numpy import float32, uint16, int32, frombuffer, dtype

COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "D": "D", "I": "I"}

def complement(nucleotide):
    """
    Complement a single nucleotide. Complements of D(eletion) and I(nsertion) are D and I, respectively.

    Args:
        nucleotide (string) : Nucleotide, must be A, C, T, G, D, or I

    Returns:
        str : Complemented nucleotide

    Raises:
        ValueError: Nucleotide must be one of A, C, T, G, D, or I
    """
    if nucleotide in COMPLEMENT_MAP:
        return COMPLEMENT_MAP[nucleotide]
    raise ValueError("Nucleotide must be one of A, C, T, G, D, or I")

def read_char(handle):
    """
    Helper function to parse character from file handle

    Args:
        handle (file): File handle

    Returns:
        char
    """
    return handle.read(1)

def read_ushort(handle):
    """
    Helper function to parse ushort from file handle

    Args:
        handle (file): File handle

    Returns:
        numpy.int16
    """
    return frombuffer(handle.read(2), dtype=uint16)[0]

def read_int(handle):
    """
    Helper function to parse int from file handle

    Args:
        handle (file): File handle

    Returns:
        numpy.int32
    """
    return struct.unpack("<i", handle.read(4))[0]

def read_float(handle):
    """
    Helper function to parse float from file handle

    Args:
        handle (file): File handle

    Returns:
        numpy.float32
    """
    return frombuffer(handle.read(4), dtype=float32)[0]

def read_byte(handle):
    """
    Helper function to parse byte from file handle

    Args:
        handle (file): File handle

    Returns:
        byte
    """
    return struct.unpack("<B", handle.read(1))[0]

def read_string(handle):
    """
    Helper function to parse string from file handle. See https://msdn.microsoft.com/en-us/library/yzxa6408(v=vs.100).aspx
    for additional details on string format.

    Args:
        handle (file): File handle

    Returns:
        string

    Raises:
        Exception: Failed to read complete string
    """
    total_length = 0
    partial_length = read_byte(handle)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << (7 * num_bytes)
        partial_length = ord(struct.unpack("c", handle.read(1))[0])
        num_bytes += 1
    total_length += partial_length << (7 * num_bytes)
    result = handle.read(total_length)
    result = result.decode("utf-8")
    if len(result) < total_length:
        raise Exception("Failed to read complete string")
    else:
        return result
