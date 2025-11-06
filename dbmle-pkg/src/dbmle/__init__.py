# src/dbmle/__init__.py

# Re-export the new core API that now lives under the dbmle package.
from dbmle.core import (
    dbmle,
    dbmle_from_ZD,
    DBMLEResult,
)

# --- Back-compat aliases (old names) ---------------------------------
# These keep existing user code working without changes.
counting_defiers_command = dbmle
counting_defiers_from_ZD = dbmle_from_ZD
CountingDefiersResult = DBMLEResult

__all__ = [
    "dbmle",
    "dbmle_from_ZD",
    "DBMLEResult",
]
