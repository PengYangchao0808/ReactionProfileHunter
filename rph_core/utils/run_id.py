from __future__ import annotations

import logging
import uuid

logger = logging.getLogger(__name__)

RUN_ID_FIELD = "run_id"


def new_run_id() -> str:
    return str(uuid.uuid4())


def is_valid_run_id(value: object) -> bool:
    if not isinstance(value, str):
        return False
    try:
        uuid.UUID(value)
    except (AttributeError, TypeError, ValueError):
        return False
    return True
