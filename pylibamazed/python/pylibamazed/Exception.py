# ============================================================================
#
# This file is part of: AMAZED
#
# Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
#
# https://www.lam.fr/
#
# This software is a computer program whose purpose is to estimate the
# spectrocopic redshift of astronomical sources (galaxy/quasar/star)
# from there 1D spectrum.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# ============================================================================
import functools
import sys
import traceback
from os.path import basename

from pylibamazed.redshift import AmzException, ErrorCode

AmazedError = AmzException  # for amazed clients compatibility reason, should be removed later


class APIException(AmzException):
    def __init__(self, errCode, message):
        filename = ""
        name = ""
        line = 0
        exc_type, exc_value, exc_traceback = sys.exc_info()
        if exc_traceback is not None:
            frame = traceback.extract_tb(exc_traceback)[-1]
            filename = basename(frame.filename)
            name = frame.name
            line = frame.lineno
        super().__init__(errCode.value,
                         message,
                         filename,
                         name,
                         line)

    @classmethod
    def fromException(cls, exception):
        return cls(ErrorCode.PYTHON_API_ERROR, str(exception))


# decorator to convert any non-amazed exception to AmzException
#  with the optional ErrorLogging
def exception_decorator(fn=None, *, logging=False):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except AmzException as e:  # will catch both AmzException and derived
                if logging:
                    e.LogError()
                raise e from None
            except Exception as e:
                amze = APIException.fromException(e)
                if logging:
                    amze.LogError()
                raise amze from e
        return wrapper
    if fn:
        return decorator(fn)
    return decorator
