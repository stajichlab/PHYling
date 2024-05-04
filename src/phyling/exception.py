"""Exception class that used in PHYling."""


class BinaryNotFoundError(OSError):
    """The binary is not installed or not found in path."""

    pass


class SeqtypeError(TypeError):
    """Invalid or ambiguous Seqtype."""

    pass


class IdenticalKeyError(KeyError):
    """Attempting to add a key that already exists to a dictionary."""

    pass


class EmptyWarning(UserWarning):
    """An empty item returned which will cause potential error in the following step."""

    pass
