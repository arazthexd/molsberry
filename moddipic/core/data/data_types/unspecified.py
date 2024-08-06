from typing import Any

from ..abstract import Data, Representation
from ..representations import UnspecifiedRep

class UnspecifiedData(Data):
    def __init__(self, init_rep: Representation | Any | None = None):
        if not isinstance(init_rep, Representation) and init_rep is not None:
            init_rep = UnspecifiedRep(init_rep)
        super().__init__(init_rep)