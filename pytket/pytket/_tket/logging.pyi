from typing import ClassVar

class level:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    critical: ClassVar[level] = ...
    debug: ClassVar[level] = ...
    err: ClassVar[level] = ...
    info: ClassVar[level] = ...
    off: ClassVar[level] = ...
    trace: ClassVar[level] = ...
    warn: ClassVar[level] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

def set_level(arg0: level) -> None: ...