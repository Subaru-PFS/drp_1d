import numpy as np
import pytest
from pylibamazed.Container import Container


# Tests in case T in a numpy array
def test_container_equality():
    # Two empty containers are equal
    assert (Container() == Container())

    array1 = np.array([1, 2, 3])
    array2 = np.array([0, 2, 3])
    # Container 1 included in container 2: not equal
    assert (Container(**{"key": array1}) != Container(**{"key": array1, "key2": array1}))

    # Containers with same data but different key not equal
    assert (Container(**{"key": array1}) != Container(**{"other key": array1}))

    # Container with same keys but different data not equal
    assert (Container(**{"key": array1}) != Container(**{"other key": array2}))

    # Containers with same data and keys but in different order are not equal
    assert (Container(**{"key": array1, "key2": array2}) == Container(**{"key2": array2, "key": array1}))


def test_container_check_append_type():
    cont = Container()
    integer = 1
    floating = 2.0

    cont.append(integer, 'key1')
    with pytest.raises(TypeError) as terr:
        cont.append(floating, 'key2')

    error_message = f"Data for container must be of the declared {type(integer)} type"
    # Attempting to append data of different type raises type error
    assert error_message in str(terr.value)
