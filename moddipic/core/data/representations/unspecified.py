from typing import Any

from ..abstract import Representation

class UnspecifiedRep(Representation):
    rep_name = "unspecified"
    def __init__(self, data: Any):
        super().__init__(data=data)