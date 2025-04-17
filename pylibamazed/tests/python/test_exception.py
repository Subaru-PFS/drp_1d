import pytest
from pylibamazed.Exception import APIException, AmzException, exception_decorator
from pylibamazed.redshift import ErrorCode


def raise_APIException():
    raise APIException(ErrorCode.PYTHON_API_ERROR, "APIException raised")


def raise_AmzException():
    raise AmzException(
        ErrorCode.INTERNAL_ERROR.value, "AmzException raised", "test_exception.py", "raise_AmzException", 11
    )


def raise_Exception():
    raise Exception("Exception raised")


@exception_decorator(logging=False)
def decorated_raise_APIException():
    raise_APIException()


@exception_decorator(logging=True)
def decorated_raise_APIException_with_logging():
    raise_APIException()


@exception_decorator(logging=False)
def decorated_raise_AmzException():
    raise_AmzException()


@exception_decorator(logging=True)
def decorated_raise_AmzException_with_logging():
    raise_AmzException()


@exception_decorator(logging=False)
def decorated_raise_Exception():
    raise_Exception()


@exception_decorator(logging=False)
def double_decorated_raise_Exception():
    decorated_raise_Exception()


@exception_decorator(logging=True)
def decorated_raise_Exception_with_logging():
    raise_Exception()


filename = "test_exception.py"


def test_APIException():
    exception = APIException
    msg = "APIException raised"
    method = "raise_APIException"
    line = 7
    short_str = "PYTHON_API_ERROR: " + msg
    long_str = short_str + f" [{filename}:{line}:{method}]"

    with pytest.raises(exception, match=msg) as exc_info:
        raise_APIException()
    e = exc_info.value
    assert e.getFileName() == ""
    assert e.getMethod() == ""
    assert e.getLine() == 0
    assert e.__str__() == short_str

    with pytest.raises(exception, match=msg) as exc_info:
        raise_APIException()
    e = exc_info.value
    e.add_traceback_info()
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_APIException()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_APIException_with_logging()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str


def test_AmzException():
    exception = AmzException
    msg = "AmzException raised"
    method = "raise_AmzException"
    line = 11
    short_str = "INTERNAL_ERROR: " + msg
    long_str = short_str + f" [{filename}:{line}:{method}]"

    with pytest.raises(exception, match=msg) as exc_info:
        raise_AmzException()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_AmzException()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_AmzException_with_logging()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str


def test_exception():
    exception = APIException
    msg = "Exception raised"
    method = "raise_Exception"
    line = 17
    short_str = "PYTHON_API_ERROR: " + "Exception: " + msg
    long_str = short_str + f" [{filename}:{line}:{method}]"

    with pytest.raises(exception, match=msg) as exc_info:
        try:
            raise_Exception()
        except Exception as ex:
            APIException.fromException(ex)
            raise APIException.fromException(ex) from None
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_Exception()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        double_decorated_raise_Exception()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(exception, match=msg) as exc_info:
        decorated_raise_Exception_with_logging()
    e = exc_info.value
    assert e.getFileName() == filename
    assert e.getMethod() == method
    assert e.getLine() == line
    assert e.__str__() == long_str

    with pytest.raises(TypeError) as exc_info:
        decorated_raise_Exception(0)  # bad call (number of arguments)
